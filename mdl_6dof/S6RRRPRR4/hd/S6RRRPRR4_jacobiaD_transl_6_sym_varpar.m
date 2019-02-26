% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:17:43
% EndTime: 2019-02-26 22:17:44
% DurationCPUTime: 0.44s
% Computational Cost: add. (718->88), mult. (599->120), div. (0->0), fcn. (470->12), ass. (0->72)
t297 = pkin(11) + qJ(5);
t293 = qJ(6) + t297;
t288 = sin(t293);
t289 = cos(t293);
t300 = qJ(2) + qJ(3);
t295 = cos(t300);
t299 = qJD(2) + qJD(3);
t347 = t295 * t299;
t294 = sin(t300);
t298 = qJD(5) + qJD(6);
t349 = t294 * t298;
t366 = t288 * t347 + t289 * t349;
t292 = cos(t297);
t281 = pkin(5) * t292 + cos(pkin(11)) * pkin(4) + pkin(3);
t356 = r_i_i_C(1) * t289;
t365 = t281 + t356;
t301 = sin(qJ(2));
t351 = pkin(2) * qJD(2);
t334 = t301 * t351;
t291 = sin(t297);
t350 = pkin(5) * qJD(5);
t336 = t291 * t350;
t296 = -pkin(10) - pkin(9) - qJ(4);
t352 = r_i_i_C(3) - t296;
t357 = pkin(5) * t291;
t364 = (t352 * t299 - t336) * t295 + (t357 + sin(pkin(11)) * pkin(4) + pkin(8) + pkin(7)) * qJD(1) - (t281 * t299 - qJD(4)) * t294 - t334;
t333 = t288 * t349;
t363 = r_i_i_C(1) * t333 + t366 * r_i_i_C(2) + qJD(4) * t295 + t294 * t336;
t304 = cos(qJ(1));
t322 = t295 * t298 - qJD(1);
t362 = t304 * t322;
t341 = qJD(1) * t295;
t321 = -t298 + t341;
t302 = sin(qJ(1));
t346 = t299 * t302;
t330 = t294 * t346;
t360 = t321 * t304 - t330;
t358 = pkin(2) * t301;
t355 = r_i_i_C(2) * t288;
t354 = r_i_i_C(2) * t289;
t353 = r_i_i_C(3) * t295;
t348 = t295 * t296;
t345 = t299 * t304;
t329 = t294 * t345;
t309 = t321 * t302 + t329;
t261 = t309 * t288 - t289 * t362;
t262 = t288 * t362 + t309 * t289;
t344 = t261 * r_i_i_C(1) + t262 * r_i_i_C(2);
t315 = t322 * t302;
t263 = t360 * t288 + t289 * t315;
t264 = t288 * t315 - t360 * t289;
t343 = -t263 * r_i_i_C(1) + t264 * r_i_i_C(2);
t340 = qJD(1) * t302;
t339 = qJD(1) * t304;
t337 = t294 * t355;
t335 = t292 * t350;
t328 = t294 * t340;
t327 = t294 * t339;
t318 = -qJD(5) + t341;
t317 = -r_i_i_C(1) * t288 - t354;
t316 = t365 * t299;
t314 = t292 * (-qJD(5) * t295 + qJD(1));
t313 = -t294 * t365 - t348;
t312 = t296 * t330 + t363 * t302 + t327 * t355 + t339 * t353;
t311 = t296 * t329 + t363 * t304 + t365 * t328 + t340 * t348;
t310 = (-r_i_i_C(3) * t294 - t295 * t365) * t299;
t303 = cos(qJ(2));
t308 = t335 + (-pkin(2) * t303 - t281 * t295 - t352 * t294 - pkin(1)) * qJD(1);
t307 = -t303 * t351 + t310;
t306 = t299 * t337 + r_i_i_C(3) * t347 + (-t296 * t299 + t317 * t298 - t336) * t295 + (qJD(4) - t316) * t294;
t274 = r_i_i_C(2) * t333;
t1 = [t264 * r_i_i_C(1) + t263 * r_i_i_C(2) - t364 * t302 + t308 * t304 (-t337 - t353 + t358) * t340 + t307 * t304 + t311 (-r_i_i_C(3) * t345 - t340 * t355) * t294 + (-r_i_i_C(3) * t340 - t304 * t316) * t295 + t311, t295 * t345 - t328 (t304 * t314 + (t318 * t302 + t329) * t291) * pkin(5) + t344, t344; -t262 * r_i_i_C(1) + t261 * r_i_i_C(2) + t308 * t302 + t364 * t304, t307 * t302 + (t313 - t358) * t339 + t312, t302 * t310 + t313 * t339 + t312, t295 * t346 + t327 (t302 * t314 + (-t318 * t304 + t330) * t291) * pkin(5) + t343, t343; 0, t306 - t334, t306, t299 * t294, t274 + (-t298 * t356 - t335) * t294 + (t317 - t357) * t347, -t366 * r_i_i_C(1) - t347 * t354 + t274;];
JaD_transl  = t1;
