% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR13
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR13_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:05
% EndTime: 2019-02-26 22:23:06
% DurationCPUTime: 0.40s
% Computational Cost: add. (362->67), mult. (867->135), div. (0->0), fcn. (1196->14), ass. (0->70)
t267 = sin(qJ(2));
t268 = sin(qJ(1));
t271 = cos(qJ(2));
t272 = cos(qJ(1));
t296 = cos(pkin(6));
t276 = t272 * t296;
t251 = t267 * t276 + t268 * t271;
t266 = sin(qJ(3));
t270 = cos(qJ(3));
t250 = t268 * t267 - t271 * t276;
t264 = cos(pkin(7));
t262 = sin(pkin(7));
t263 = sin(pkin(6));
t287 = t263 * t272;
t278 = t262 * t287;
t274 = t250 * t264 + t278;
t232 = -t251 * t270 + t274 * t266;
t243 = -t250 * t262 + t264 * t287;
t261 = pkin(13) + qJ(5);
t259 = sin(t261);
t260 = cos(t261);
t222 = t232 * t260 + t243 * t259;
t265 = sin(qJ(6));
t300 = t222 * t265;
t269 = cos(qJ(6));
t299 = t222 * t269;
t220 = t232 * t259 - t243 * t260;
t295 = t251 * t266;
t293 = t259 * t262;
t292 = t260 * t262;
t291 = t260 * t265;
t290 = t260 * t269;
t289 = t262 * t263;
t288 = t263 * t268;
t286 = t264 * t266;
t285 = t264 * t270;
t284 = t266 * t267;
t283 = t266 * t271;
t282 = t267 * t270;
t281 = t270 * t271;
t280 = t267 * t289;
t279 = t262 * t288;
t277 = t268 * t296;
t275 = t296 * t262;
t252 = -t272 * t267 - t271 * t277;
t273 = -t252 * t262 + t264 * t288;
t253 = -t267 * t277 + t272 * t271;
t249 = t296 * t264 - t271 * t289;
t248 = (-t264 * t284 + t281) * t263;
t247 = (t264 * t282 + t283) * t263;
t242 = t266 * t275 + (t264 * t283 + t282) * t263;
t241 = -t270 * t275 + (-t264 * t281 + t284) * t263;
t239 = t252 * t270 - t253 * t286;
t238 = t252 * t266 + t253 * t285;
t237 = -t250 * t270 - t251 * t286;
t236 = -t250 * t266 + t251 * t285;
t235 = t248 * t260 + t259 * t280;
t234 = t253 * t270 + (t252 * t264 + t279) * t266;
t233 = -t252 * t285 + t253 * t266 - t270 * t279;
t231 = -t274 * t270 - t295;
t229 = t250 * t285 + t270 * t278 + t295;
t228 = t242 * t260 + t249 * t259;
t227 = -t242 * t259 + t249 * t260;
t226 = t239 * t260 + t253 * t293;
t225 = t237 * t260 + t251 * t293;
t224 = t234 * t260 + t273 * t259;
t223 = t234 * t259 - t273 * t260;
t219 = t224 * t269 + t233 * t265;
t218 = -t224 * t265 + t233 * t269;
t1 = [t231 * t265 + t299, t226 * t269 + t238 * t265, -t233 * t290 + t234 * t265, 0, -t223 * t269, t218; t219, t225 * t269 + t236 * t265, -t229 * t290 - t232 * t265, 0, t220 * t269, t229 * t269 + t300; 0, t235 * t269 + t247 * t265, -t241 * t290 + t242 * t265, 0, t227 * t269, -t228 * t265 + t241 * t269; t231 * t269 - t300, -t226 * t265 + t238 * t269, t233 * t291 + t234 * t269, 0, t223 * t265, -t219; t218, -t225 * t265 + t236 * t269, t229 * t291 - t232 * t269, 0, -t220 * t265, -t229 * t265 + t299; 0, -t235 * t265 + t247 * t269, t241 * t291 + t242 * t269, 0, -t227 * t265, -t228 * t269 - t241 * t265; t220, t239 * t259 - t253 * t292, -t233 * t259, 0, t224, 0; t223, t237 * t259 - t251 * t292, -t229 * t259, 0, -t222, 0; 0, t248 * t259 - t260 * t280, -t241 * t259, 0, t228, 0;];
JR_rot  = t1;
