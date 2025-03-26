from flask import Flask, render_template
from pages.pagina1 import pagina1_bp
from pages.pagina2 import pagina2_bp

app = Flask(__name__)
app.secret_key = 'tu_clave_secreta_aqui'


@app.route('/')
def main():
    return render_template('main.html')

# Registrar blueprints
app.register_blueprint(pagina2_bp)
app.register_blueprint(pagina1_bp)

if __name__ == '__main__':
    app.run(debug=True)
