# Personal website built with Hugo
(theme: [Hugo-classic](https://github.com/goodroot/hugo-classic))

### Hugo and Github Pages
To use a custom domain we must create a file `static/CNAME` that only contains the CNAME of the custom domain (which, in this case, is www.johannesnauta.com). See: https://gohugo.io/hosting-and-deployment/hosting-on-github/
Then, just push to the `main` branch, and Github will set in motion what needs to happen. The details on the Github Action are in `.github/workflows/gh-pages.yml`.

### Adding content
Most content is in the `/content/` folder. Images, or other static items, can be put into `/static/`. Please find other examples, such as `_index.md` how to add an image.

##### Adding specific/new content
...

#### Overriding CSS
The theme contains css in `themes/hugo-classic/`
