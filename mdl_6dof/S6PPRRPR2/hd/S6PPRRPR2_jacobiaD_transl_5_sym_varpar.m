% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRPR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRPR2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRPR2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:58
% EndTime: 2019-02-26 19:40:59
% DurationCPUTime: 0.22s
% Computational Cost: add. (338->54), mult. (1082->105), div. (0->0), fcn. (1228->12), ass. (0->52)
t325 = sin(pkin(12));
t328 = sin(pkin(6));
t334 = sin(qJ(3));
t336 = cos(qJ(3));
t329 = cos(pkin(12));
t331 = cos(pkin(7));
t352 = t329 * t331;
t327 = sin(pkin(7));
t332 = cos(pkin(6));
t355 = t327 * t332;
t312 = (t325 * t336 + t334 * t352) * t328 + t334 * t355;
t330 = cos(pkin(11));
t326 = sin(pkin(11));
t357 = t326 * t332;
t321 = -t325 * t357 + t330 * t329;
t320 = -t330 * t325 - t329 * t357;
t356 = t327 * t328;
t340 = t320 * t331 + t326 * t356;
t308 = t321 * t336 + t340 * t334;
t351 = t330 * t332;
t319 = t325 * t351 + t326 * t329;
t318 = -t326 * t325 + t329 * t351;
t348 = t330 * t356;
t341 = -t318 * t331 + t348;
t364 = -t319 * t336 + t341 * t334;
t363 = pkin(4) - r_i_i_C(2);
t362 = -pkin(9) - r_i_i_C(1);
t361 = r_i_i_C(3) + qJ(5);
t360 = t319 * t334;
t354 = t328 * t329;
t353 = t328 * t331;
t350 = qJD(3) * t334;
t349 = qJD(3) * t336;
t346 = t327 * t349;
t345 = t331 * t349;
t313 = -t318 * t327 - t330 * t353;
t333 = sin(qJ(4));
t335 = cos(qJ(4));
t344 = t313 * t333 - t335 * t364;
t314 = -t320 * t327 + t326 * t353;
t343 = t308 * t335 + t314 * t333;
t317 = -t327 * t354 + t332 * t331;
t342 = t312 * t335 + t317 * t333;
t338 = qJD(3) * (t361 * t333 + t363 * t335 + pkin(3));
t337 = qJD(5) * t333 + (-t363 * t333 + t361 * t335) * qJD(4);
t309 = t328 * t325 * t350 - t332 * t346 - t345 * t354;
t303 = -t326 * t328 * t346 - t320 * t345 + t321 * t350;
t301 = -t318 * t345 + (t336 * t348 + t360) * qJD(3);
t299 = t342 * qJD(4) - t309 * t333;
t297 = t343 * qJD(4) - t303 * t333;
t295 = t344 * qJD(4) - t301 * t333;
t1 = [0, 0, t362 * t303 + t337 * (-t321 * t334 + t340 * t336) - t308 * t338, t343 * qJD(5) + t361 * (-t303 * t335 + (-t308 * t333 + t314 * t335) * qJD(4)) - t363 * t297, t297, 0; 0, 0, t362 * t301 + t337 * (-t341 * t336 - t360) + t364 * t338, t344 * qJD(5) + t361 * (-t301 * t335 + (t313 * t335 + t333 * t364) * qJD(4)) - t363 * t295, t295, 0; 0, 0, t362 * t309 + t337 * (t336 * t355 + (-t325 * t334 + t336 * t352) * t328) - t312 * t338, t342 * qJD(5) + t361 * (-t309 * t335 + (-t312 * t333 + t317 * t335) * qJD(4)) - t363 * t299, t299, 0;];
JaD_transl  = t1;
