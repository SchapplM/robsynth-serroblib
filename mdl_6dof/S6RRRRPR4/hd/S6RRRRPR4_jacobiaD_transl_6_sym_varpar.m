% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:32:18
% EndTime: 2019-02-26 22:32:18
% DurationCPUTime: 0.42s
% Computational Cost: add. (740->88), mult. (636->118), div. (0->0), fcn. (492->12), ass. (0->71)
t304 = qJ(4) + pkin(11);
t297 = qJ(6) + t304;
t292 = sin(t297);
t293 = cos(t297);
t305 = qJ(2) + qJ(3);
t299 = cos(t305);
t303 = qJD(2) + qJD(3);
t352 = t299 * t303;
t298 = sin(t305);
t302 = qJD(4) + qJD(6);
t354 = t298 * t302;
t369 = t292 * t352 + t293 * t354;
t285 = pkin(5) * cos(t304) + cos(qJ(4)) * pkin(4);
t283 = pkin(3) + t285;
t360 = r_i_i_C(1) * t293;
t368 = t283 + t360;
t284 = sin(qJ(4)) * pkin(4) + pkin(5) * sin(t304);
t280 = t284 * qJD(4);
t307 = sin(qJ(2));
t355 = pkin(2) * qJD(2);
t340 = t307 * t355;
t301 = -pkin(10) - qJ(5) - pkin(9);
t356 = r_i_i_C(3) - t301;
t367 = (t356 * t303 - t280) * t299 + (t284 + pkin(8) + pkin(7)) * qJD(1) - (t283 * t303 - qJD(5)) * t298 - t340;
t339 = t292 * t354;
t366 = r_i_i_C(1) * t339 + t369 * r_i_i_C(2) + qJD(5) * t299;
t311 = cos(qJ(1));
t328 = t299 * t302 - qJD(1);
t365 = t311 * t328;
t345 = qJD(1) * t299;
t327 = -t302 + t345;
t308 = sin(qJ(1));
t350 = t303 * t308;
t336 = t298 * t350;
t363 = t327 * t311 - t336;
t361 = pkin(2) * t307;
t359 = r_i_i_C(2) * t292;
t358 = r_i_i_C(2) * t293;
t357 = r_i_i_C(3) * t299;
t353 = t299 * t301;
t351 = t303 * t298;
t349 = t303 * t311;
t335 = t298 * t349;
t317 = t327 * t308 + t335;
t262 = t317 * t292 - t293 * t365;
t263 = t292 * t365 + t317 * t293;
t348 = t262 * r_i_i_C(1) + t263 * r_i_i_C(2);
t322 = t328 * t308;
t264 = t363 * t292 + t293 * t322;
t265 = t292 * t322 - t363 * t293;
t347 = -t264 * r_i_i_C(1) + t265 * r_i_i_C(2);
t344 = qJD(1) * t308;
t343 = qJD(1) * t311;
t341 = t298 * t359;
t334 = t298 * t344;
t333 = t298 * t343;
t325 = -r_i_i_C(1) * t292 - t358;
t324 = t284 * t345 - t280;
t323 = t368 * t303;
t321 = -t341 - t357;
t320 = t301 * t336 + t366 * t308 + t333 * t359 + t343 * t357;
t319 = t301 * t335 + t366 * t311 + t368 * t334 + t344 * t353;
t281 = t285 * qJD(4);
t318 = qJD(1) * t285 - t281 * t299 + t284 * t351;
t310 = cos(qJ(2));
t316 = t281 + (-pkin(2) * t310 - t283 * t299 - t356 * t298 - pkin(1)) * qJD(1);
t315 = t280 * t298 + (-r_i_i_C(3) * t298 - t299 * t368) * t303;
t314 = -t310 * t355 + t315;
t313 = r_i_i_C(3) * t352 + t303 * t341 + (-t301 * t303 + t325 * t302 - t280) * t299 + (qJD(5) - t323) * t298;
t276 = r_i_i_C(2) * t339;
t1 = [t265 * r_i_i_C(1) + t264 * r_i_i_C(2) - t367 * t308 + t316 * t311 (t321 + t361) * t344 + t314 * t311 + t319, t315 * t311 + t321 * t344 + t319, t324 * t308 + t318 * t311 + t348, t299 * t349 - t334, t348; -t263 * r_i_i_C(1) + t262 * r_i_i_C(2) + t316 * t308 + t367 * t311 (-t298 * t368 - t353 - t361) * t343 + t314 * t308 + t320 (-t301 * t343 - t308 * t323) * t299 + ((-r_i_i_C(3) * t303 + t280) * t308 - t368 * t343) * t298 + t320, t318 * t308 - t324 * t311 + t347, t299 * t350 + t333, t347; 0, t313 - t340, t313, t276 + (-t302 * t360 - t281) * t298 + (-t284 + t325) * t352, t351, -t369 * r_i_i_C(1) - t352 * t358 + t276;];
JaD_transl  = t1;
