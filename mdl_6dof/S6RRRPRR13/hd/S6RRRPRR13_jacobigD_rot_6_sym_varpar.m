% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR13_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:06
% EndTime: 2019-02-26 22:23:06
% DurationCPUTime: 0.18s
% Computational Cost: add. (112->50), mult. (332->104), div. (0->0), fcn. (355->12), ass. (0->48)
t314 = sin(pkin(7));
t319 = sin(qJ(2));
t348 = t314 * t319;
t315 = sin(pkin(6));
t320 = sin(qJ(1));
t347 = t315 * t320;
t323 = cos(qJ(1));
t346 = t315 * t323;
t316 = cos(pkin(7));
t318 = sin(qJ(3));
t345 = t316 * t318;
t344 = t318 * t319;
t322 = cos(qJ(2));
t343 = t318 * t322;
t321 = cos(qJ(3));
t342 = t319 * t321;
t341 = t320 * t319;
t340 = t320 * t322;
t339 = t321 * t322;
t338 = t323 * t319;
t337 = t323 * t322;
t336 = qJD(1) * t315;
t313 = pkin(13) + qJ(5);
t311 = sin(t313);
t335 = qJD(3) * t311;
t334 = qJD(3) * t314;
t333 = t320 * t336;
t332 = t323 * t336;
t331 = t314 * t333;
t330 = t314 * t332;
t317 = cos(pkin(6));
t307 = t317 * t337 - t341;
t329 = t307 * t316 - t314 * t346;
t309 = -t317 * t340 - t338;
t328 = t309 * t316 + t314 * t347;
t327 = t316 * t343 + t342;
t308 = t317 * t338 + t340;
t326 = t317 * t341 - t337;
t325 = t308 * t321 + t329 * t318;
t324 = t328 * t318 - t321 * t326;
t312 = cos(t313);
t306 = -t326 * qJD(1) + t307 * qJD(2);
t305 = t309 * qJD(1) - t308 * qJD(2);
t304 = -t308 * qJD(1) + t309 * qJD(2);
t303 = -t307 * qJD(1) + t326 * qJD(2);
t302 = -t305 * t314 + t316 * t333;
t301 = -t303 * t314 + t316 * t332;
t1 = [0, t332, t301, 0, t304 * t318 + (-t303 * t316 - t330) * t321 + t324 * qJD(3) (t303 * t345 + t304 * t321 + t318 * t330) * t311 - t301 * t312 + (t324 * t312 + (-t309 * t314 + t316 * t347) * t311) * qJD(5) + (t318 * t326 + t328 * t321) * t335; 0, t333, t302, 0, t306 * t318 + (-t305 * t316 - t331) * t321 + t325 * qJD(3) (t305 * t345 + t306 * t321 + t318 * t331) * t311 - t302 * t312 + (t325 * t312 + (-t307 * t314 - t316 * t346) * t311) * qJD(5) + (-t308 * t318 + t329 * t321) * t335; 0, 0, t315 * qJD(2) * t348, 0, t317 * t318 * t334 + (t327 * qJD(3) + (t316 * t342 + t343) * qJD(2)) * t315 (t321 * t311 * t334 + (t312 * t314 * t318 + t311 * t316) * qJD(5)) * t317 + ((-t314 * t322 * t311 + t327 * t312) * qJD(5) + (t316 * t339 - t344) * t335 + ((-t316 * t344 + t339) * t311 - t312 * t348) * qJD(2)) * t315;];
JgD_rot  = t1;
