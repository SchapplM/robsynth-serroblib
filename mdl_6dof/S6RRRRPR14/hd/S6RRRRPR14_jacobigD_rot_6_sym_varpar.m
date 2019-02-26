% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPR14_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:12
% EndTime: 2019-02-26 22:38:13
% DurationCPUTime: 0.22s
% Computational Cost: add. (100->49), mult. (332->104), div. (0->0), fcn. (355->12), ass. (0->46)
t304 = sin(pkin(7));
t312 = cos(qJ(4));
t339 = t304 * t312;
t305 = sin(pkin(6));
t311 = sin(qJ(1));
t338 = t305 * t311;
t315 = cos(qJ(1));
t337 = t305 * t315;
t306 = cos(pkin(7));
t309 = sin(qJ(3));
t336 = t306 * t309;
t310 = sin(qJ(2));
t335 = t309 * t310;
t314 = cos(qJ(2));
t334 = t309 * t314;
t313 = cos(qJ(3));
t333 = t310 * t313;
t332 = t311 * t310;
t331 = t311 * t314;
t330 = t313 * t314;
t329 = t315 * t310;
t328 = t315 * t314;
t327 = qJD(1) * t305;
t308 = sin(qJ(4));
t326 = qJD(3) * t308;
t325 = t311 * t327;
t324 = t315 * t327;
t323 = t304 * t325;
t322 = t304 * t324;
t307 = cos(pkin(6));
t300 = t307 * t328 - t332;
t321 = t300 * t306 - t304 * t337;
t302 = -t307 * t331 - t329;
t320 = t302 * t306 + t304 * t338;
t319 = t306 * t334 + t333;
t301 = t307 * t329 + t331;
t318 = t307 * t332 - t328;
t317 = t301 * t313 + t321 * t309;
t316 = t320 * t309 - t313 * t318;
t299 = -t318 * qJD(1) + t300 * qJD(2);
t298 = t302 * qJD(1) - t301 * qJD(2);
t297 = -t301 * qJD(1) + t302 * qJD(2);
t296 = -t300 * qJD(1) + t318 * qJD(2);
t295 = -t298 * t304 + t306 * t325;
t294 = -t296 * t304 + t306 * t324;
t1 = [0, t324, t294, t297 * t309 + (-t296 * t306 - t322) * t313 + t316 * qJD(3), 0 (t296 * t336 + t297 * t313 + t309 * t322) * t308 - t294 * t312 + (t316 * t312 + (-t302 * t304 + t306 * t338) * t308) * qJD(4) + (t309 * t318 + t320 * t313) * t326; 0, t325, t295, t299 * t309 + (-t298 * t306 - t323) * t313 + t317 * qJD(3), 0 (t298 * t336 + t299 * t313 + t309 * t323) * t308 - t295 * t312 + (t317 * t312 + (-t300 * t304 - t306 * t337) * t308) * qJD(4) + (-t301 * t309 + t321 * t313) * t326; 0, 0, t305 * qJD(2) * t310 * t304, t307 * t304 * qJD(3) * t309 + (t319 * qJD(3) + (t306 * t333 + t334) * qJD(2)) * t305, 0 (t304 * t313 * t326 + (t306 * t308 + t309 * t339) * qJD(4)) * t307 + ((-t304 * t314 * t308 + t319 * t312) * qJD(4) + (t306 * t330 - t335) * t326 + ((-t306 * t335 + t330) * t308 - t310 * t339) * qJD(2)) * t305;];
JgD_rot  = t1;
