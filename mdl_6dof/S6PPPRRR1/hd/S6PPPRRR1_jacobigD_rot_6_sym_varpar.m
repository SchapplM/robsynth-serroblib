% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPPRRR1_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobigD_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:45
% EndTime: 2019-02-26 19:38:45
% DurationCPUTime: 0.12s
% Computational Cost: add. (117->35), mult. (361->80), div. (0->0), fcn. (460->16), ass. (0->45)
t297 = sin(pkin(12));
t306 = cos(pkin(6));
t325 = t297 * t306;
t299 = sin(pkin(7));
t300 = sin(pkin(6));
t324 = t299 * t300;
t323 = t299 * t306;
t305 = cos(pkin(7));
t322 = t300 * t305;
t302 = cos(pkin(13));
t321 = t302 * t305;
t303 = cos(pkin(12));
t320 = t303 * t306;
t307 = sin(qJ(5));
t319 = qJD(4) * t307;
t296 = sin(pkin(13));
t292 = t296 * t320 + t297 * t302;
t295 = sin(pkin(14));
t301 = cos(pkin(14));
t291 = -t296 * t297 + t302 * t320;
t315 = t291 * t305 - t303 * t324;
t282 = -t292 * t295 + t301 * t315;
t288 = -t291 * t299 - t303 * t322;
t298 = sin(pkin(8));
t304 = cos(pkin(8));
t318 = t282 * t304 + t288 * t298;
t294 = -t296 * t325 + t302 * t303;
t293 = -t296 * t303 - t302 * t325;
t314 = t293 * t305 + t297 * t324;
t284 = -t294 * t295 + t301 * t314;
t289 = -t293 * t299 + t297 * t322;
t317 = t284 * t304 + t289 * t298;
t286 = t301 * t323 + (-t295 * t296 + t301 * t321) * t300;
t290 = -t302 * t324 + t305 * t306;
t316 = t286 * t304 + t290 * t298;
t283 = t292 * t301 + t295 * t315;
t308 = sin(qJ(4));
t310 = cos(qJ(4));
t313 = t283 * t310 + t308 * t318;
t285 = t294 * t301 + t295 * t314;
t312 = t285 * t310 + t308 * t317;
t287 = t300 * t296 * t301 + (t300 * t321 + t323) * t295;
t311 = t287 * t310 + t308 * t316;
t309 = cos(qJ(5));
t1 = [0, 0, 0, 0, t312 * qJD(4) (t312 * t309 + (-t284 * t298 + t289 * t304) * t307) * qJD(5) + (-t285 * t308 + t310 * t317) * t319; 0, 0, 0, 0, t313 * qJD(4) (t313 * t309 + (-t282 * t298 + t288 * t304) * t307) * qJD(5) + (-t283 * t308 + t310 * t318) * t319; 0, 0, 0, 0, t311 * qJD(4) (t311 * t309 + (-t286 * t298 + t290 * t304) * t307) * qJD(5) + (-t287 * t308 + t310 * t316) * t319;];
JgD_rot  = t1;
