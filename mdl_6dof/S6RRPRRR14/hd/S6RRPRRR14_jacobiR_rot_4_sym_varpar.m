% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR14_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiR_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:21
% EndTime: 2018-12-10 18:38:22
% DurationCPUTime: 0.28s
% Computational Cost: add. (738->74), mult. (731->114), div. (0->0), fcn. (746->26), ass. (0->74)
t279 = sin(qJ(2));
t280 = sin(qJ(1));
t283 = cos(qJ(1));
t296 = pkin(6) + qJ(2);
t289 = cos(t296) / 0.2e1;
t297 = pkin(6) - qJ(2);
t293 = cos(t297);
t284 = t293 / 0.2e1 + t289;
t239 = t280 * t279 - t283 * t284;
t287 = sin(t296) / 0.2e1;
t291 = sin(t297);
t253 = t287 - t291 / 0.2e1;
t282 = cos(qJ(2));
t240 = t253 * t283 + t280 * t282;
t268 = pkin(7) + pkin(14);
t257 = sin(t268) / 0.2e1;
t269 = pkin(7) - pkin(14);
t266 = sin(t269);
t246 = t257 + t266 / 0.2e1;
t258 = cos(t269) / 0.2e1;
t267 = cos(t268);
t248 = t258 + t267 / 0.2e1;
t270 = sin(pkin(14));
t273 = sin(pkin(6));
t298 = t273 * t283;
t223 = t239 * t248 + t240 * t270 + t246 * t298;
t247 = t257 - t266 / 0.2e1;
t249 = t258 - t267 / 0.2e1;
t274 = cos(pkin(14));
t224 = t239 * t247 - t240 * t274 + t249 * t298;
t272 = sin(pkin(7));
t276 = cos(pkin(7));
t235 = -t239 * t272 + t276 * t298;
t294 = pkin(8) + qJ(4);
t286 = sin(t294) / 0.2e1;
t295 = pkin(8) - qJ(4);
t290 = sin(t295);
t250 = t286 + t290 / 0.2e1;
t288 = cos(t294) / 0.2e1;
t292 = cos(t295);
t255 = t292 / 0.2e1 + t288;
t278 = sin(qJ(4));
t305 = t223 * t255 - t224 * t278 + t235 * t250;
t251 = t286 - t290 / 0.2e1;
t254 = t288 - t292 / 0.2e1;
t281 = cos(qJ(4));
t304 = t223 * t251 + t224 * t281 - t235 * t254;
t303 = t250 * t272;
t302 = t254 * t272;
t256 = t289 - t293 / 0.2e1;
t301 = t256 * t272;
t275 = cos(pkin(8));
t300 = t272 * t275;
t299 = t273 * t280;
t244 = t280 * t253 - t282 * t283;
t242 = -t283 * t279 - t280 * t284;
t225 = t242 * t248 + t244 * t270 + t246 * t299;
t226 = t242 * t247 - t244 * t274 + t249 * t299;
t237 = -t242 * t272 + t276 * t299;
t285 = t225 * t251 + t226 * t281 - t237 * t254;
t277 = cos(pkin(6));
t271 = sin(pkin(8));
t252 = t287 + t291 / 0.2e1;
t238 = -t252 * t272 + t276 * t277;
t234 = t247 * t256 + t252 * t274;
t233 = t248 * t256 - t252 * t270;
t232 = t247 * t252 + t249 * t277 - t256 * t274;
t231 = t246 * t277 + t248 * t252 + t256 * t270;
t230 = t242 * t274 + t244 * t247;
t229 = -t242 * t270 + t244 * t248;
t228 = -t239 * t274 - t240 * t247;
t227 = t239 * t270 - t240 * t248;
t220 = t225 * t255 - t226 * t278 + t237 * t250;
t1 = [t304, t229 * t251 + t230 * t281 + t244 * t302, 0, t220, 0, 0; t285, t227 * t251 + t228 * t281 - t240 * t302, 0, -t305, 0, 0; 0, t233 * t251 + t234 * t281 + t254 * t301, 0, t231 * t255 - t232 * t278 + t238 * t250, 0, 0; t305, t229 * t255 - t230 * t278 - t244 * t303, 0, -t285, 0, 0; t220, t227 * t255 - t228 * t278 + t240 * t303, 0, t304, 0, 0; 0, t233 * t255 - t234 * t278 - t250 * t301, 0, -t231 * t251 - t232 * t281 + t238 * t254, 0, 0; -t223 * t271 + t235 * t275, -t229 * t271 - t244 * t300, 0, 0, 0, 0; -t225 * t271 + t237 * t275, -t227 * t271 + t240 * t300, 0, 0, 0, 0; 0, -t233 * t271 - t256 * t300, 0, 0, 0, 0;];
JR_rot  = t1;
