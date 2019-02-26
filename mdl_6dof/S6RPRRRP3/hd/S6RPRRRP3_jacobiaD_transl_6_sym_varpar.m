% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:09:04
% EndTime: 2019-02-26 21:09:04
% DurationCPUTime: 0.35s
% Computational Cost: add. (739->71), mult. (742->105), div. (0->0), fcn. (622->10), ass. (0->58)
t302 = qJ(4) + qJ(5);
t299 = cos(t302);
t346 = r_i_i_C(3) + qJ(6);
t359 = t346 * t299;
t300 = qJD(4) + qJD(5);
t298 = sin(t302);
t303 = sin(qJ(4));
t345 = pkin(4) * qJD(4);
t347 = r_i_i_C(2) + pkin(9) + pkin(8);
t310 = t347 * qJD(3) + qJD(6) * t298 - t303 * t345;
t349 = pkin(5) + r_i_i_C(1);
t333 = t349 * t298;
t358 = (t333 - t359) * t300 - t310;
t305 = cos(qJ(4));
t295 = pkin(4) * t305 + pkin(3);
t306 = cos(qJ(3));
t304 = sin(qJ(3));
t337 = qJD(3) * t304;
t357 = -t295 * t337 + t310 * t306;
t316 = -t346 * t298 - t349 * t299;
t313 = -t295 + t316;
t355 = t313 * t304 + t347 * t306;
t354 = (qJD(1) * t306 - qJD(4)) * t303;
t343 = t298 * t306;
t350 = t299 * t337 + t300 * t343;
t348 = pkin(4) * t303;
t344 = t298 * t300;
t342 = t299 * t300;
t341 = t299 * t306;
t340 = t300 * t304;
t301 = qJ(1) + pkin(10);
t296 = sin(t301);
t339 = qJD(1) * t296;
t297 = cos(t301);
t338 = qJD(1) * t297;
t336 = qJD(3) * t306;
t335 = qJD(6) * t299;
t334 = t305 * t345;
t332 = t297 * t342;
t330 = t299 * t340;
t329 = t304 * t335 + t336 * t359;
t328 = pkin(7) + t348;
t327 = t299 * t339;
t326 = t298 * t337;
t320 = t334 - t335;
t319 = t296 * t298 + t297 * t341;
t318 = -t295 * t306 - t347 * t304 - pkin(2);
t317 = t296 * t342 + t298 * t338;
t276 = t297 * t326 - t306 * t332 - t296 * t344 + (t296 * t343 + t297 * t299) * qJD(1);
t277 = t350 * t297 + t306 * t327 - t317;
t315 = t319 * qJD(6) + t349 * t276 - t346 * t277;
t278 = -t296 * t326 - t297 * t344 + t317 * t306 - t327;
t279 = t319 * qJD(1) - t350 * t296 - t332;
t314 = -(-t296 * t341 + t297 * t298) * qJD(6) + t346 * t279 - t349 * t278;
t312 = t303 * t337 + (-qJD(4) * t306 + qJD(1)) * t305;
t311 = qJD(3) * t313;
t308 = t358 * t304 + t306 * t311;
t1 = [t320 * t297 - t349 * t279 - t346 * t278 - t357 * t296 + (-cos(qJ(1)) * pkin(1) - t328 * t296 + t318 * t297) * qJD(1), 0, t308 * t297 - t355 * t339 (t296 * t354 + t312 * t297) * pkin(4) + t315, t315, -t276; t320 * t296 - t349 * t277 - t346 * t276 + t357 * t297 + (-sin(qJ(1)) * pkin(1) + t328 * t297 + t318 * t296) * qJD(1), 0, t308 * t296 + t355 * t338 (t312 * t296 - t297 * t354) * pkin(4) + t314, t314, t278; 0, 0, t304 * t311 - t358 * t306 (-t333 - t348) * t336 + (t316 * t300 - t334) * t304 + t329, -t349 * t330 + (-t349 * t336 - t346 * t340) * t298 + t329, t298 * t336 + t330;];
JaD_transl  = t1;
