% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR12_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobiR_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:10
% EndTime: 2019-02-26 21:21:11
% DurationCPUTime: 0.40s
% Computational Cost: add. (399->55), mult. (1173->112), div. (0->0), fcn. (1599->16), ass. (0->63)
t263 = cos(pkin(6));
t260 = cos(pkin(14));
t271 = cos(qJ(1));
t278 = t271 * t260;
t256 = sin(pkin(14));
t267 = sin(qJ(1));
t281 = t267 * t256;
t250 = -t263 * t278 + t281;
t279 = t271 * t256;
t280 = t267 * t260;
t251 = t263 * t279 + t280;
t266 = sin(qJ(3));
t270 = cos(qJ(3));
t258 = sin(pkin(7));
t259 = sin(pkin(6));
t285 = t259 * t271;
t277 = t258 * t285;
t262 = cos(pkin(7));
t282 = t262 * t266;
t240 = t250 * t282 - t251 * t270 + t266 * t277;
t265 = sin(qJ(4));
t269 = cos(qJ(4));
t239 = (t250 * t262 + t277) * t270 + t251 * t266;
t246 = -t250 * t258 + t262 * t285;
t257 = sin(pkin(8));
t261 = cos(pkin(8));
t275 = t239 * t261 + t246 * t257;
t224 = t240 * t269 + t275 * t265;
t231 = t239 * t257 - t246 * t261;
t264 = sin(qJ(5));
t268 = cos(qJ(5));
t299 = t224 * t264 + t231 * t268;
t298 = t224 * t268 - t231 * t264;
t222 = t240 * t265 - t275 * t269;
t252 = -t263 * t280 - t279;
t286 = t259 * t267;
t248 = -t252 * t258 + t262 * t286;
t291 = t248 * t257;
t289 = t257 * t264;
t288 = t257 * t268;
t287 = t258 * t263;
t284 = t261 * t265;
t283 = t261 * t269;
t244 = t270 * t287 + (t260 * t262 * t270 - t256 * t266) * t259;
t249 = -t259 * t260 * t258 + t263 * t262;
t274 = t244 * t261 + t249 * t257;
t272 = t252 * t262 + t258 * t286;
t253 = -t263 * t281 + t278;
t245 = t266 * t287 + (t256 * t270 + t260 * t282) * t259;
t242 = t253 * t270 + t272 * t266;
t241 = -t253 * t266 + t272 * t270;
t236 = -t244 * t257 + t249 * t261;
t234 = t244 * t269 - t245 * t284;
t233 = -t241 * t257 + t248 * t261;
t230 = t245 * t269 + t274 * t265;
t229 = -t245 * t265 + t274 * t269;
t228 = t241 * t269 - t242 * t284;
t227 = -t239 * t269 + t240 * t284;
t226 = t242 * t269 + (t241 * t261 + t291) * t265;
t225 = -t241 * t283 + t242 * t265 - t269 * t291;
t221 = t226 * t268 + t233 * t264;
t220 = -t226 * t264 + t233 * t268;
t1 = [t298, 0, t228 * t268 + t242 * t289, -t225 * t268, t220, 0; t221, 0, t227 * t268 - t240 * t289, t222 * t268, t299, 0; 0, 0, t234 * t268 + t245 * t289, t229 * t268, -t230 * t264 + t236 * t268, 0; -t299, 0, -t228 * t264 + t242 * t288, t225 * t264, -t221, 0; t220, 0, -t227 * t264 - t240 * t288, -t222 * t264, t298, 0; 0, 0, -t234 * t264 + t245 * t288, -t229 * t264, -t230 * t268 - t236 * t264, 0; t222, 0, t241 * t265 + t242 * t283, t226, 0, 0; t225, 0, -t239 * t265 - t240 * t283, -t224, 0, 0; 0, 0, t244 * t265 + t245 * t283, t230, 0, 0;];
JR_rot  = t1;
