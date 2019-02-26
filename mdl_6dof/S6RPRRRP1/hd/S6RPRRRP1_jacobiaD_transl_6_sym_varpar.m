% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP1
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
% Datum: 2019-02-26 21:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:07:52
% EndTime: 2019-02-26 21:07:53
% DurationCPUTime: 0.40s
% Computational Cost: add. (672->73), mult. (706->100), div. (0->0), fcn. (567->10), ass. (0->57)
t313 = qJ(1) + pkin(10);
t309 = cos(t313);
t317 = cos(qJ(5));
t308 = sin(t313);
t354 = qJD(1) * t308;
t342 = t317 * t354;
t315 = sin(qJ(5));
t352 = qJD(5) * t315;
t370 = t309 * t352 + t342;
t363 = -r_i_i_C(1) - pkin(5);
t361 = r_i_i_C(3) + qJ(6);
t369 = t361 * t315;
t314 = qJ(3) + qJ(4);
t311 = cos(t314);
t364 = pkin(9) + r_i_i_C(2);
t346 = t364 * t311;
t350 = qJD(6) * t315;
t351 = qJD(5) * t317;
t368 = t361 * t351 + t350;
t312 = qJD(3) + qJD(4);
t316 = sin(qJ(3));
t360 = pkin(3) * qJD(3);
t348 = t316 * t360;
t310 = sin(t314);
t358 = t310 * t312;
t366 = -pkin(4) * t358 + (t364 * t312 + t350) * t311 - t348;
t362 = pkin(3) * t316;
t359 = t308 * t315;
t357 = t311 * t312;
t356 = t311 * t317;
t355 = t312 * t315;
t353 = qJD(1) * t309;
t349 = t317 * qJD(6);
t347 = t364 * t310;
t345 = t310 * t355;
t344 = t317 * t358;
t340 = t311 * t352;
t338 = t308 * t352;
t336 = t309 * t351;
t335 = -t363 * t310 * t338 + t353 * t346;
t330 = t309 * t356 + t359;
t318 = cos(qJ(3));
t328 = -t318 * pkin(3) - pkin(4) * t311 - pkin(2) - t347;
t327 = (-t363 * t370 + (pkin(4) + t369) * t354) * t310;
t326 = t363 * t317 - t369;
t325 = t308 * t351 + t315 * t353;
t324 = -pkin(4) + t326;
t323 = t324 * t310;
t322 = t368 * t311 + t312 * t323 + t363 * t340 + t364 * t357;
t321 = -t368 * t310 + (t324 * t311 - t347) * t312;
t320 = -t318 * t360 + t321;
t319 = -pkin(8) - pkin(7);
t280 = t330 * qJD(1) - t308 * t344 - t311 * t338 - t336;
t279 = -t308 * t345 + t325 * t311 - t370;
t278 = t311 * t342 + (t340 + t344) * t309 - t325;
t277 = t309 * t345 - t311 * t336 - t338 + (t309 * t317 + t311 * t359) * qJD(1);
t1 = [-t309 * t349 + t363 * t280 - t361 * t279 - t366 * t308 + (-cos(qJ(1)) * pkin(1) + t308 * t319 + t328 * t309) * qJD(1), 0 (-t346 + t362) * t354 + t320 * t309 + t327, t321 * t309 - t346 * t354 + t327, t330 * qJD(6) - t363 * t277 - t361 * t278, -t277; -t308 * t349 + t363 * t278 - t361 * t277 + t366 * t309 + (-sin(qJ(1)) * pkin(1) - t309 * t319 + t328 * t308) * qJD(1), 0 (t323 - t362) * t353 + t320 * t308 + t335, t321 * t308 + t323 * t353 + t335 -(-t308 * t356 + t309 * t315) * qJD(6) + t361 * t280 + t363 * t279, t279; 0, 0, t322 - t348, t322 (t363 * t315 + t361 * t317) * t357 + (t326 * qJD(5) + t349) * t310, t310 * t351 + t311 * t355;];
JaD_transl  = t1;
