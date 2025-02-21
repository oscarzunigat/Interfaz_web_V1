from flask import Flask, render_template
from pages.pagina1 import pagina1_bp

app = Flask(__name__)


@app.route('/')
def main():
    return render_template('main.html')

@app.route('/pagina2')
def pagina2():
    return render_template('pagina2.html')

################################ pagina 1 ################################

app.register_blueprint(pagina1_bp)

@app.route('/pagina1')
def pagina1():
    return render_template('pagina1.html')

##########################################################################
if __name__ == '__main__':
    app.run(debug=True)