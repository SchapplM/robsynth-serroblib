% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR12_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:37:00
% EndTime: 2019-02-26 22:37:00
% DurationCPUTime: 0.40s
% Computational Cost: add. (362->67), mult. (867->135), div. (0->0), fcn. (1196->14), ass. (0->70)
t268 = sin(qJ(2));
t269 = sin(qJ(1));
t272 = cos(qJ(2));
t273 = cos(qJ(1));
t297 = cos(pkin(6));
t277 = t273 * t297;
t252 = t268 * t277 + t269 * t272;
t267 = sin(qJ(3));
t271 = cos(qJ(3));
t251 = t268 * t269 - t272 * t277;
t265 = cos(pkin(7));
t263 = sin(pkin(7));
t264 = sin(pkin(6));
t288 = t264 * t273;
t279 = t263 * t288;
t275 = t251 * t265 + t279;
t233 = -t252 * t271 + t275 * t267;
t244 = -t251 * t263 + t265 * t288;
t262 = qJ(4) + pkin(13);
t260 = sin(t262);
t261 = cos(t262);
t223 = t233 * t261 + t244 * t260;
t266 = sin(qJ(6));
t301 = t223 * t266;
t270 = cos(qJ(6));
t300 = t223 * t270;
t221 = t233 * t260 - t244 * t261;
t296 = t252 * t267;
t294 = t260 * t263;
t293 = t261 * t263;
t292 = t261 * t266;
t291 = t261 * t270;
t290 = t263 * t264;
t289 = t264 * t269;
t287 = t265 * t267;
t286 = t265 * t271;
t285 = t267 * t268;
t284 = t267 * t272;
t283 = t268 * t271;
t282 = t271 * t272;
t281 = t268 * t290;
t280 = t263 * t289;
t278 = t269 * t297;
t276 = t297 * t263;
t253 = -t273 * t268 - t272 * t278;
t274 = -t253 * t263 + t265 * t289;
t254 = -t268 * t278 + t272 * t273;
t250 = t297 * t265 - t272 * t290;
t249 = (-t265 * t285 + t282) * t264;
t248 = (t265 * t283 + t284) * t264;
t243 = t267 * t276 + (t265 * t284 + t283) * t264;
t242 = -t271 * t276 + (-t265 * t282 + t285) * t264;
t240 = t253 * t271 - t254 * t287;
t239 = t253 * t267 + t254 * t286;
t238 = -t251 * t271 - t252 * t287;
t237 = -t251 * t267 + t252 * t286;
t236 = t249 * t261 + t260 * t281;
t235 = t254 * t271 + (t253 * t265 + t280) * t267;
t234 = -t253 * t286 + t254 * t267 - t271 * t280;
t232 = -t275 * t271 - t296;
t230 = t251 * t286 + t271 * t279 + t296;
t229 = t243 * t261 + t250 * t260;
t228 = -t243 * t260 + t250 * t261;
t227 = t240 * t261 + t254 * t294;
t226 = t238 * t261 + t252 * t294;
t225 = t235 * t261 + t274 * t260;
t224 = t235 * t260 - t274 * t261;
t220 = t225 * t270 + t234 * t266;
t219 = -t225 * t266 + t234 * t270;
t1 = [t232 * t266 + t300, t227 * t270 + t239 * t266, -t234 * t291 + t235 * t266, -t224 * t270, 0, t219; t220, t226 * t270 + t237 * t266, -t230 * t291 - t233 * t266, t221 * t270, 0, t230 * t270 + t301; 0, t236 * t270 + t248 * t266, -t242 * t291 + t243 * t266, t228 * t270, 0, -t229 * t266 + t242 * t270; t232 * t270 - t301, -t227 * t266 + t239 * t270, t234 * t292 + t235 * t270, t224 * t266, 0, -t220; t219, -t226 * t266 + t237 * t270, t230 * t292 - t233 * t270, -t221 * t266, 0, -t230 * t266 + t300; 0, -t236 * t266 + t248 * t270, t242 * t292 + t243 * t270, -t228 * t266, 0, -t229 * t270 - t242 * t266; t221, t240 * t260 - t254 * t293, -t234 * t260, t225, 0, 0; t224, t238 * t260 - t252 * t293, -t230 * t260, -t223, 0, 0; 0, t249 * t260 - t261 * t281, -t242 * t260, t229, 0, 0;];
JR_rot  = t1;
