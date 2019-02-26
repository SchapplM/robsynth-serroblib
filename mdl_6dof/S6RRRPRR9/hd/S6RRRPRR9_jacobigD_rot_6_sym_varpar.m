% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR9
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
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR9_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:38
% EndTime: 2019-02-26 22:20:39
% DurationCPUTime: 0.23s
% Computational Cost: add. (143->58), mult. (464->118), div. (0->0), fcn. (512->14), ass. (0->48)
t327 = sin(pkin(13));
t330 = cos(pkin(13));
t334 = sin(qJ(3));
t338 = cos(qJ(3));
t343 = t327 * t334 - t330 * t338;
t356 = t343 * qJD(3);
t328 = sin(pkin(7));
t335 = sin(qJ(2));
t355 = t328 * t335;
t329 = sin(pkin(6));
t336 = sin(qJ(1));
t354 = t329 * t336;
t340 = cos(qJ(1));
t353 = t329 * t340;
t352 = t335 * t340;
t351 = t336 * t335;
t339 = cos(qJ(2));
t350 = t336 * t339;
t349 = t339 * t340;
t348 = qJD(1) * t336;
t347 = qJD(1) * t340;
t346 = t329 * t348;
t345 = t329 * t347;
t344 = t327 * t338 + t330 * t334;
t332 = cos(pkin(6));
t321 = t332 * t349 - t351;
t323 = -t332 * t350 - t352;
t322 = t332 * t352 + t350;
t342 = t332 * t351 - t349;
t320 = t344 * qJD(3);
t337 = cos(qJ(5));
t333 = sin(qJ(5));
t331 = cos(pkin(7));
t318 = t344 * t331;
t317 = t343 * t331;
t316 = t344 * t328;
t315 = t343 * t328;
t314 = t331 * t356;
t313 = t331 * t320;
t312 = t328 * t356;
t311 = t328 * t320;
t310 = -t342 * qJD(1) + t321 * qJD(2);
t309 = t323 * qJD(1) - t322 * qJD(2);
t308 = -t322 * qJD(1) + t323 * qJD(2);
t307 = -t321 * qJD(1) + t342 * qJD(2);
t306 = -t309 * t328 + t331 * t346;
t305 = -t307 * t328 + t331 * t345;
t1 = [0, t345, t305, 0, t307 * t317 + t308 * t344 + t323 * t313 + t342 * t356 + (t311 * t336 + t315 * t347) * t329 (t307 * t318 - t308 * t343 - t323 * t314 + t342 * t320 + (-t312 * t336 + t316 * t347) * t329) * t333 - t305 * t337 + ((t316 * t354 + t318 * t323 + t342 * t343) * t337 + (-t323 * t328 + t331 * t354) * t333) * qJD(5); 0, t346, t306, 0, t309 * t317 + t310 * t344 + t321 * t313 - t322 * t356 + (-t311 * t340 + t315 * t348) * t329 (t309 * t318 - t310 * t343 - t321 * t314 - t322 * t320 + (t312 * t340 + t316 * t348) * t329) * t333 - t306 * t337 + ((-t316 * t353 + t321 * t318 - t322 * t343) * t337 + (-t321 * t328 - t331 * t353) * t333) * qJD(5); 0, 0, t329 * qJD(2) * t355, 0, t311 * t332 + (t313 * t339 - t356 * t335 + (-t317 * t335 + t339 * t344) * qJD(2)) * t329 (-t312 * t333 + (t316 * t337 + t331 * t333) * qJD(5)) * t332 + ((-t314 * t339 - t320 * t335) * t333 + ((t318 * t339 - t335 * t343) * t337 - t328 * t339 * t333) * qJD(5) + ((-t318 * t335 - t339 * t343) * t333 - t337 * t355) * qJD(2)) * t329;];
JgD_rot  = t1;
