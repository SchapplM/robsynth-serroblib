% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR9
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
% Datum: 2019-02-22 12:08
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR9_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:08:26
% EndTime: 2019-02-22 12:08:26
% DurationCPUTime: 0.36s
% Computational Cost: add. (461->63), mult. (1310->135), div. (0->0), fcn. (1807->16), ass. (0->65)
t272 = sin(pkin(7));
t271 = sin(pkin(13));
t277 = sin(qJ(3));
t300 = cos(pkin(13));
t302 = cos(qJ(3));
t287 = t271 * t302 + t277 * t300;
t258 = t287 * t272;
t274 = cos(pkin(7));
t260 = t287 * t274;
t278 = sin(qJ(2));
t279 = sin(qJ(1));
t282 = cos(qJ(2));
t283 = cos(qJ(1));
t301 = cos(pkin(6));
t290 = t283 * t301;
t262 = t279 * t278 - t282 * t290;
t263 = t278 * t290 + t279 * t282;
t286 = -t277 * t271 + t300 * t302;
t273 = sin(pkin(6));
t295 = t273 * t283;
t239 = t258 * t295 + t262 * t260 - t263 * t286;
t254 = -t262 * t272 + t274 * t295;
t276 = sin(qJ(5));
t281 = cos(qJ(5));
t228 = t239 * t281 + t254 * t276;
t275 = sin(qJ(6));
t280 = cos(qJ(6));
t257 = t286 * t272;
t259 = t286 * t274;
t288 = -t257 * t295 - t262 * t259 - t263 * t287;
t306 = t228 * t275 - t280 * t288;
t305 = t228 * t280 + t275 * t288;
t226 = t239 * t276 - t254 * t281;
t299 = t272 * t273;
t298 = t272 * t276;
t297 = t272 * t281;
t296 = t273 * t279;
t294 = t275 * t281;
t293 = t280 * t281;
t292 = t278 * t299;
t291 = t279 * t301;
t264 = -t283 * t278 - t282 * t291;
t289 = -t264 * t272 + t274 * t296;
t265 = -t278 * t291 + t283 * t282;
t285 = t258 * t296 + t264 * t260 + t265 * t286;
t284 = t301 * t258 + (t260 * t282 + t278 * t286) * t273;
t261 = t274 * t301 - t282 * t299;
t252 = (-t260 * t278 + t282 * t286) * t273;
t251 = (t259 * t278 + t282 * t287) * t273;
t250 = t252 * t281 + t276 * t292;
t249 = -t265 * t260 + t264 * t286;
t248 = t265 * t259 + t264 * t287;
t247 = -t263 * t260 - t262 * t286;
t246 = t263 * t259 - t262 * t287;
t244 = t301 * t257 + (t259 * t282 - t278 * t287) * t273;
t241 = t257 * t296 + t264 * t259 - t265 * t287;
t234 = t249 * t281 + t265 * t298;
t233 = t247 * t281 + t263 * t298;
t232 = t261 * t276 + t281 * t284;
t231 = t261 * t281 - t276 * t284;
t230 = t276 * t289 + t281 * t285;
t229 = t276 * t285 - t281 * t289;
t225 = t230 * t280 - t241 * t275;
t224 = -t230 * t275 - t241 * t280;
t1 = [t305, t234 * t280 + t248 * t275, t241 * t293 + t275 * t285, 0, -t229 * t280, t224; t225, t233 * t280 + t246 * t275, -t239 * t275 + t288 * t293, 0, t226 * t280, t306; 0, t250 * t280 + t251 * t275, t244 * t293 + t275 * t284, 0, t231 * t280, -t232 * t275 - t244 * t280; -t306, -t234 * t275 + t248 * t280, -t241 * t294 + t280 * t285, 0, t229 * t275, -t225; t224, -t233 * t275 + t246 * t280, -t239 * t280 - t288 * t294, 0, -t226 * t275, t305; 0, -t250 * t275 + t251 * t280, -t244 * t294 + t280 * t284, 0, -t231 * t275, -t232 * t280 + t244 * t275; t226, t249 * t276 - t265 * t297, t241 * t276, 0, t230, 0; t229, t247 * t276 - t263 * t297, t288 * t276, 0, -t228, 0; 0, t252 * t276 - t281 * t292, t244 * t276, 0, t232, 0;];
JR_rot  = t1;
