% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:40:43
% EndTime: 2019-02-26 22:40:43
% DurationCPUTime: 0.37s
% Computational Cost: add. (695->94), mult. (678->129), div. (0->0), fcn. (522->10), ass. (0->75)
t301 = qJ(4) + qJ(5);
t293 = sin(t301);
t308 = cos(qJ(1));
t299 = qJD(4) + qJD(5);
t302 = qJ(2) + qJ(3);
t296 = cos(t302);
t343 = qJD(1) * t296;
t326 = -t299 + t343;
t294 = sin(t302);
t300 = qJD(2) + qJD(3);
t305 = sin(qJ(1));
t348 = t300 * t305;
t335 = t294 * t348;
t364 = t326 * t308 - t335;
t369 = t364 * t293;
t295 = cos(t301);
t306 = cos(qJ(4));
t285 = t306 * pkin(4) + pkin(5) * t295;
t283 = pkin(3) + t285;
t368 = r_i_i_C(1) * t295 + t283;
t303 = sin(qJ(4));
t354 = pkin(4) * qJD(4);
t360 = pkin(5) * t293;
t274 = -t299 * t360 - t303 * t354;
t284 = t303 * pkin(4) + t360;
t304 = sin(qJ(2));
t355 = pkin(2) * qJD(2);
t338 = t304 * t355;
t298 = -qJ(6) - pkin(10) - pkin(9);
t356 = r_i_i_C(3) - t298;
t367 = (t356 * t300 + t274) * t296 + (t284 + pkin(8) + pkin(7)) * qJD(1) - (t283 * t300 - qJD(6)) * t294 - t338;
t350 = t296 * t300;
t336 = t293 * t350;
t353 = t294 * t299;
t337 = t293 * t353;
t352 = t295 * t299;
t358 = r_i_i_C(2) * t294;
t366 = r_i_i_C(1) * t337 + r_i_i_C(2) * t336 + qJD(6) * t296 + t352 * t358;
t362 = -pkin(5) - r_i_i_C(1);
t361 = pkin(2) * t304;
t357 = r_i_i_C(3) * t296;
t351 = t296 * t298;
t349 = t300 * t294;
t347 = t300 * t308;
t334 = t294 * t347;
t314 = t326 * t305 + t334;
t327 = t296 * t299 - qJD(1);
t321 = t295 * t327;
t262 = t314 * t293 - t308 * t321;
t263 = t327 * t308 * t293 + t314 * t295;
t346 = t262 * r_i_i_C(1) + t263 * r_i_i_C(2);
t320 = t327 * t305;
t264 = t295 * t320 + t369;
t265 = t293 * t320 - t295 * t364;
t345 = -t264 * r_i_i_C(1) + t265 * r_i_i_C(2);
t342 = qJD(1) * t305;
t341 = qJD(1) * t308;
t339 = t293 * t358;
t333 = t294 * t342;
t332 = t294 * t341;
t324 = -r_i_i_C(1) * t293 - r_i_i_C(2) * t295;
t323 = t284 * t343 + t274;
t322 = t368 * t300;
t319 = -t339 - t357;
t318 = t293 * r_i_i_C(2) * t332 + t298 * t335 + t366 * t305 + t341 * t357;
t317 = t298 * t334 + t366 * t308 + t368 * t333 + t342 * t351;
t275 = pkin(5) * t352 + t306 * t354;
t316 = qJD(1) * t285 - t275 * t296 + t284 * t349;
t307 = cos(qJ(2));
t313 = t275 + (-t307 * pkin(2) - t283 * t296 - t356 * t294 - pkin(1)) * qJD(1);
t312 = -t274 * t294 + (-r_i_i_C(3) * t294 - t296 * t368) * t300;
t311 = -t307 * t355 + t312;
t310 = r_i_i_C(3) * t350 + t300 * t339 + (-t298 * t300 + t324 * t299 + t274) * t296 + (qJD(6) - t322) * t294;
t280 = r_i_i_C(2) * t337;
t1 = [t265 * r_i_i_C(1) + t264 * r_i_i_C(2) - t367 * t305 + t313 * t308 (t319 + t361) * t342 + t311 * t308 + t317, t312 * t308 + t319 * t342 + t317, t323 * t305 + t316 * t308 + t346, t262 * pkin(5) + t346, t296 * t347 - t333; -t263 * r_i_i_C(1) + t262 * r_i_i_C(2) + t313 * t305 + t367 * t308 (-t294 * t368 - t351 - t361) * t341 + t311 * t305 + t318 (-t298 * t341 - t305 * t322) * t296 + ((-r_i_i_C(3) * t300 - t274) * t305 - t368 * t341) * t294 + t318, t316 * t305 - t323 * t308 + t345 (-t305 * t321 - t369) * pkin(5) + t345, t296 * t348 + t332; 0, t310 - t338, t310, t280 + (-r_i_i_C(1) * t352 - t275) * t294 + (-t284 + t324) * t350, t280 + t362 * t336 + (-r_i_i_C(2) * t350 + t362 * t353) * t295, t349;];
JaD_transl  = t1;
