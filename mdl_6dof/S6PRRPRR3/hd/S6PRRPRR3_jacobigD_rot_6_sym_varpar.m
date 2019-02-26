% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRR3_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:15
% EndTime: 2019-02-26 20:05:15
% DurationCPUTime: 0.20s
% Computational Cost: add. (109->49), mult. (360->103), div. (0->0), fcn. (404->14), ass. (0->40)
t290 = sin(pkin(13));
t294 = cos(pkin(13));
t299 = sin(qJ(3));
t302 = cos(qJ(3));
t306 = t299 * t290 - t302 * t294;
t314 = t306 * qJD(3);
t291 = sin(pkin(12));
t293 = sin(pkin(6));
t313 = t291 * t293;
t292 = sin(pkin(7));
t301 = cos(qJ(5));
t312 = t292 * t301;
t295 = cos(pkin(12));
t311 = t293 * t295;
t296 = cos(pkin(7));
t310 = t293 * t296;
t297 = cos(pkin(6));
t300 = sin(qJ(2));
t309 = t297 * t300;
t303 = cos(qJ(2));
t308 = t297 * t303;
t307 = t302 * t290 + t299 * t294;
t282 = -t291 * t300 + t295 * t308;
t283 = t291 * t303 + t295 * t309;
t284 = -t291 * t308 - t295 * t300;
t305 = t291 * t309 - t295 * t303;
t287 = t307 * qJD(3);
t298 = sin(qJ(5));
t281 = t305 * qJD(2);
t280 = t284 * qJD(2);
t279 = t283 * qJD(2);
t278 = t282 * qJD(2);
t277 = t307 * t296;
t276 = t306 * t296;
t275 = t307 * t292;
t274 = t296 * t314;
t273 = t296 * t287;
t272 = t292 * t314;
t271 = t292 * t287;
t1 = [0, 0, -t281 * t292, 0, t271 * t313 + t284 * t273 + t281 * t276 + t280 * t307 + t305 * t314 (-t272 * t313 - t284 * t274 + t281 * t277 - t280 * t306 + t287 * t305) * t298 + t281 * t312 + ((t275 * t313 + t284 * t277 + t305 * t306) * t301 + (-t284 * t292 + t291 * t310) * t298) * qJD(5); 0, 0, t279 * t292, 0, -t271 * t311 + t282 * t273 - t279 * t276 + t278 * t307 - t283 * t314 (t272 * t311 - t282 * t274 - t279 * t277 - t278 * t306 - t283 * t287) * t298 - t279 * t312 + ((-t275 * t311 + t282 * t277 - t283 * t306) * t301 + (-t282 * t292 - t295 * t310) * t298) * qJD(5); 0, 0, t293 * qJD(2) * t300 * t292, 0, t297 * t271 + (t273 * t303 - t314 * t300 + (-t276 * t300 + t303 * t307) * qJD(2)) * t293 (-t272 * t298 + (t275 * t301 + t296 * t298) * qJD(5)) * t297 + ((-t274 * t303 - t287 * t300) * t298 + ((t277 * t303 - t300 * t306) * t301 - t292 * t303 * t298) * qJD(5) + ((-t277 * t300 - t303 * t306) * t298 - t300 * t312) * qJD(2)) * t293;];
JgD_rot  = t1;
