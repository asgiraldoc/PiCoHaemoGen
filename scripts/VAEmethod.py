# Imported libraries
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow.keras.layers import *
from tensorflow.keras.models import Model
from tensorflow.keras.metrics import categorical_crossentropy
from tensorflow.keras import backend as K

# Disabling eager execution for TF 1.x compatibility
tf.compat.v1.disable_eager_execution()

def store_data(file_name, data_obj):
    """Stores data from a dictionary into a file."""
    with open(file_name, "w") as output_file:
        for i in range(data_obj["shape"][0]):
            print(data_obj["name"][i], data_obj["group"][i], end="", file=output_file)
            for j in range(data_obj["shape"][1]):
                print(" ", end=" ", file=output_file)
                print(",".join([str(val) for val in list(data_obj["dna_code"][i][j])]), end="", file=output_file)
            print("", file=output_file)


def retrieve_data(file_name):
    """Loads data from a file into a dictionary."""
    data_obj = {
        "data_set": file_name.replace(".txt", ""),
        "name": [],
        "group": [],
        "dna_code": []
    }
    with open(file_name, "r") as input_file:
        for line in input_file:
            values = line.rstrip("\n").split(" ")
            data_obj["name"].append(values[0])
            data_obj["group"].append(values[1])
            data_obj["dna_code"].append([val.split(",") for val in values[2:]])

    data_obj["name"] = np.array(data_obj["name"], str)
    data_obj["group"] = np.array(data_obj["group"], str)
    data_obj["dna_code"] = np.array(data_obj["dna_code"], float)
    data_obj["shape"] = data_obj["dna_code"].shape
    return data_obj


def plot_z_distribution(z, loc, title="none", legend=True):
    """Plots the Z distribution."""
    x = z[:, 0]
    y = z[:, 1]

    for idx, lo in enumerate(np.unique(loc)):
        xx = []
        yy = []
        for i in range(z.shape[0]):
            if lo == loc[i]:
                xx.append(x[i])
                yy.append(y[i])
        plt.scatter(xx, yy, label=lo, color="#1f77b4")

    if title != "none": plt.title(title)
    if legend == True:
        plt.legend(bbox_to_anchor=(1, 0, 0.5, 1), loc="upper left",)


def plot_mu_sg(mu, sg, loc, sample=100, alpha=0.5, title="none", legend=True):
    """Plots the mean and standard deviation."""
    for idx, lo in enumerate(np.unique(loc)):
        xx = []
        yy = []
        for i in range(len(loc)):
            if lo == loc[i]:
                x, y = (mu[i] + np.random.normal(0, 1, size=(sample, 2)) * sg[i]).T
                xx += list(x)
                yy += list(y)
        plt.scatter(xx, yy, label=lo, alpha=alpha, color="#1f77b4")

    plt.scatter(mu[:, 0], mu[:, 1], s=30, facecolors='none', edgecolors='black', linewidths=1)
    if title != "none": plt.title(title)
    if legend == True:
        plt.legend(bbox_to_anchor=(1, 0, 0.5, 1), loc="upper left",)


def mk_model(
    original_dim,    # number of snps
    cat,             # number of categories of one-hot-encoding
    latent_dim=2,
    
    # encoder
    en_dim=[100,100,100],
    en_drop=[0.5,0.5,0.5],
    
    # decoder
    de_dim=[100,100,100],  # number of neurons for each layer
    de_drop=[0.5,0.5,0.5], # rate of dropout for each layer
    
    act = "elu" # activation function for each layer
):

  def act_fn(fn,tensor):
    if fn == "leakyrelu": return LeakyReLU()(tensor)
    else: return Activation(fn)(tensor)
  
  half_cat = int(cat/2)
  ########################################################################
  # INPUT
  ########################################################################  
  x_in = Input(shape=(original_dim,cat),name="x_in")
  x_in_em = Dense(half_cat,use_bias=False,name="x_in_em")(x_in)
    
  ########################################################################
  # ENCODER :: Q(z|X)
  ########################################################################
  en = Flatten()(x_in_em)
  en = BatchNormalization(scale=False,center=False)(en)
  
  for i in range(len(en_dim)):
    en = Dense(en_dim[i])(en)
    en = Dropout(en_drop[i])(en)
    en = act_fn(act,en)
    en = BatchNormalization(scale=False,center=False)(en)
    
  ########################################################################
  # Z (Latent space)
  ########################################################################  
  Z_mu = Dense(latent_dim)(en)
  Z_log_sigma_sq = Dense(latent_dim)(en)
  
  Z_sigma = Lambda(lambda x: K.exp(0.5*x))(Z_log_sigma_sq)
  
  Z = Lambda(lambda x: x[0]+x[1]*K.random_normal(K.shape(x[0])))([Z_mu,Z_sigma])
  
  ########################################################################
  # DECODER :: P(X|z)
  ########################################################################
  de = Z
  for i in range(len(de_dim)):
    de = Dense(de_dim[i])(de)
    de = Dropout(de_drop[i])(de)
    de = act_fn(act,de)
    de = BatchNormalization(scale=False,center=False)(de)
  
  de = Dense(original_dim*half_cat)(de)
  ########################################################################
  # OUTPUT
  ########################################################################  
  x_out_em = Reshape((-1,half_cat))(de)
  x_out = Dense(cat,activation="softmax")(x_out_em)
  ########################################################################
  
  def vae_loss(kl_weight=0.5):
    def loss(x_true, x_pred):
      # mask out missing data!
      mask = K.sum(x_in,axis=-1)

      # sigma (or standard deviation), keeping close to 1
      # mu (or mean), keeping close to 0
      kl_loss = K.sum(K.square(Z_mu) + K.square(Z_sigma) - Z_log_sigma_sq - 
                      1.0, axis=-1)

      # reconstruction (categorical crossentropy)
      recon = K.sum(categorical_crossentropy(x_in,x_out) * mask, axis=-1)
      
      return K.mean(recon + kl_loss * kl_weight)
    return loss
    
  def acc(x_true,x_pred):
    mask = K.sum(x_in,axis=-1,keepdims=True)
    acc = K.sum(K.square(x_in-x_out),axis=-1,keepdims=True)
    return K.mean(1.0 - K.sqrt(K.sum(acc*mask,axis=1)/K.sum(mask,axis=1)))
    
  vae0 = Model([x_in],[x_out],name="vae0")
  vae0.compile(optimizer='adam', loss=vae_loss(0.1), metrics=[acc])
  
  vae1 = Model([x_in],[x_out],name="vae1")
  vae1.compile(optimizer='adam', loss=vae_loss(0.5), metrics=[acc])
  
  enc = Model([x_in],[Z_mu,Z_sigma],name="enc")
  return vae0,vae1,enc

def do_it(data):
  plt.rcParams['figure.figsize'] = [10, 10]
  plt.style.use('seaborn-colorblind')

  def gen(batch_size):
    while True:
      idx = np.random.randint(0,data["shape"][0],size=batch_size)
      tmp = data["dna_code"][idx]
      yield tmp,tmp

  K.clear_session()

  vae0,vae1,enc = mk_model(data["shape"][1],data["shape"][2])
  loss_history = []
  acc_history = []
  r = 4
  for i in range(r):
    f = 1/(r-i)
    batch_size = int(data["shape"][0] * f + 0.5)
    steps = int(data["shape"][0]/batch_size + 0.5)
    epochs = int(1000 * f + 0.5)

    for vae in (vae0,vae1):
      print("-")
      his = vae.fit(
          gen(batch_size),
          steps_per_epoch=steps,
          epochs=epochs,
          verbose=False
      )
      loss_history += list(his.history['loss'])
      acc_history += list(his.history['acc'])

    if i == r-1:
      plt.subplot(2, 2, 1)
      plt.plot(np.arange(len(loss_history)),loss_history)
      plt.title("Loss function training")
      plt.xlabel("Epochs")
      plt.ylabel("Loss")
      plt.subplot(2, 2, 2)
      plt.plot(np.arange(len(acc_history)),acc_history)
      plt.title("Accuracy function training")
      plt.ylabel("Accuracy")
      plt.xlabel("Epochs")

      vae_mu,vae_sg = enc.predict(data["dna_code"])
      plt.subplot(2, 2, 3)
      plot_z_distribution(vae_mu,data["group"],legend=False)
      plt.title("Z distribution")
      plt.subplot(2, 2, 4)
      plot_mu_sg(vae_mu,vae_sg,data["group"],sample=100,legend=False)
      plt.title(r'$\mu$ and $\sigma$ distribution')
      plt.savefig(data["data_set"]+".png")
      
  return vae_mu,vae_sg, data 


def VAE_model(input_file):
  mu,sg,data = do_it(retrieve_data(input_file))
  return mu, data
