% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JgD_rot [3x7]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S7RRRRRRR1_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobigD_rot_6_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobigD_rot_6_sym_varpar: qJD has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobigD_rot_6_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:22
% EndTime: 2019-02-26 22:54:22
% DurationCPUTime: 0.18s
% Computational Cost: add. (81->44), mult. (250->87), div. (0->0), fcn. (261->10), ass. (0->43)
t279 = cos(qJ(3));
t311 = -qJD(4) * t279 + qJD(2);
t275 = sin(qJ(2));
t276 = sin(qJ(1));
t280 = cos(qJ(2));
t290 = qJD(1) * t280 + qJD(3);
t281 = cos(qJ(1));
t297 = qJD(2) * t281;
t310 = t275 * t297 + t290 * t276;
t272 = sin(qJ(5));
t273 = sin(qJ(4));
t309 = t272 * t273;
t278 = cos(qJ(4));
t308 = t272 * t278;
t307 = t273 * t275;
t306 = t276 * t279;
t305 = t276 * t280;
t274 = sin(qJ(3));
t304 = t281 * t274;
t303 = t281 * t279;
t302 = qJD(1) * t276;
t301 = qJD(1) * t281;
t300 = qJD(2) * t275;
t299 = qJD(2) * t279;
t298 = qJD(2) * t280;
t296 = qJD(3) * t274;
t295 = qJD(3) * t279;
t294 = qJD(4) * t275;
t291 = -qJD(3) * t280 - qJD(1);
t289 = -qJD(4) + t299;
t286 = t291 * t281;
t288 = t274 * t286 - t310 * t279 + t281 * t294;
t287 = t290 * t303 + (t291 * t274 - t275 * t299 + t294) * t276;
t285 = -t275 * t301 - t276 * t298;
t284 = t275 * t302 - t280 * t297;
t270 = t279 * t305 + t304;
t283 = -qJD(4) * t270 - t285;
t271 = -t276 * t274 + t280 * t303;
t282 = qJD(4) * t271 + t284;
t277 = cos(qJ(5));
t268 = t291 * t306 + (t276 * t300 - t290 * t281) * t274;
t266 = t310 * t274 + t279 * t286;
t1 = [0, t301, t284, t266, t288 * t273 + t282 * t278, -t266 * t277 + t288 * t308 - t282 * t309 + ((t271 * t278 + t281 * t307) * t277 + (-t280 * t304 - t306) * t272) * qJD(5), 0; 0, t302, t285, t268, t287 * t273 - t283 * t278, -t268 * t277 + t287 * t308 + t283 * t309 + ((t270 * t278 + t276 * t307) * t277 + (-t274 * t305 + t303) * t272) * qJD(5), 0; 0, 0, -t300, -t274 * t298 - t275 * t295, -t311 * t278 * t275 + (-t275 * t296 + t289 * t280) * t273 (t289 * t308 + (qJD(2) * t274 - t273 * qJD(5)) * t277) * t280 + ((t311 * t273 - t278 * t296) * t272 + t277 * t295 + (t278 * t279 * t277 - t274 * t272) * qJD(5)) * t275, 0;];
JgD_rot  = t1;
