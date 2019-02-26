% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP11_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:39
% EndTime: 2019-02-26 22:45:40
% DurationCPUTime: 0.42s
% Computational Cost: add. (293->66), mult. (867->135), div. (0->0), fcn. (1196->14), ass. (0->69)
t257 = sin(qJ(2));
t258 = sin(qJ(1));
t262 = cos(qJ(2));
t263 = cos(qJ(1));
t287 = cos(pkin(6));
t267 = t263 * t287;
t243 = t257 * t267 + t258 * t262;
t256 = sin(qJ(3));
t261 = cos(qJ(3));
t242 = t258 * t257 - t262 * t267;
t253 = cos(pkin(7));
t251 = sin(pkin(7));
t252 = sin(pkin(6));
t280 = t252 * t263;
t269 = t251 * t280;
t265 = t242 * t253 + t269;
t224 = -t243 * t261 + t265 * t256;
t235 = -t242 * t251 + t253 * t280;
t255 = sin(qJ(4));
t260 = cos(qJ(4));
t214 = t224 * t260 + t235 * t255;
t254 = sin(qJ(5));
t291 = t214 * t254;
t259 = cos(qJ(5));
t290 = t214 * t259;
t212 = t224 * t255 - t235 * t260;
t286 = t243 * t256;
t284 = t251 * t252;
t283 = t251 * t255;
t282 = t251 * t260;
t281 = t252 * t258;
t279 = t253 * t256;
t278 = t253 * t261;
t277 = t254 * t260;
t276 = t256 * t257;
t275 = t256 * t262;
t274 = t257 * t261;
t273 = t259 * t260;
t272 = t261 * t262;
t271 = t257 * t284;
t270 = t251 * t281;
t268 = t258 * t287;
t266 = t287 * t251;
t244 = -t263 * t257 - t262 * t268;
t264 = -t244 * t251 + t253 * t281;
t245 = -t257 * t268 + t263 * t262;
t241 = t287 * t253 - t262 * t284;
t240 = (-t253 * t276 + t272) * t252;
t239 = (t253 * t274 + t275) * t252;
t234 = t256 * t266 + (t253 * t275 + t274) * t252;
t233 = -t261 * t266 + (-t253 * t272 + t276) * t252;
t231 = t240 * t260 + t255 * t271;
t230 = t244 * t261 - t245 * t279;
t229 = t244 * t256 + t245 * t278;
t228 = -t242 * t261 - t243 * t279;
t227 = -t242 * t256 + t243 * t278;
t226 = t245 * t261 + (t244 * t253 + t270) * t256;
t225 = -t244 * t278 + t245 * t256 - t261 * t270;
t223 = -t265 * t261 - t286;
t221 = t242 * t278 + t261 * t269 + t286;
t220 = t234 * t260 + t241 * t255;
t219 = -t234 * t255 + t241 * t260;
t218 = t230 * t260 + t245 * t283;
t217 = t228 * t260 + t243 * t283;
t216 = t226 * t260 + t264 * t255;
t215 = t226 * t255 - t264 * t260;
t211 = t216 * t259 + t225 * t254;
t210 = -t216 * t254 + t225 * t259;
t1 = [t223 * t254 + t290, t218 * t259 + t229 * t254, -t225 * t273 + t226 * t254, -t215 * t259, t210, 0; t211, t217 * t259 + t227 * t254, -t221 * t273 - t224 * t254, t212 * t259, t221 * t259 + t291, 0; 0, t231 * t259 + t239 * t254, -t233 * t273 + t234 * t254, t219 * t259, -t220 * t254 + t233 * t259, 0; t223 * t259 - t291, -t218 * t254 + t229 * t259, t225 * t277 + t226 * t259, t215 * t254, -t211, 0; t210, -t217 * t254 + t227 * t259, t221 * t277 - t224 * t259, -t212 * t254, -t221 * t254 + t290, 0; 0, -t231 * t254 + t239 * t259, t233 * t277 + t234 * t259, -t219 * t254, -t220 * t259 - t233 * t254, 0; t212, t230 * t255 - t245 * t282, -t225 * t255, t216, 0, 0; t215, t228 * t255 - t243 * t282, -t221 * t255, -t214, 0, 0; 0, t240 * t255 - t260 * t271, -t233 * t255, t220, 0, 0;];
JR_rot  = t1;
