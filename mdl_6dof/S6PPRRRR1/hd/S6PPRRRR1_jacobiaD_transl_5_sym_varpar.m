% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:44
% EndTime: 2019-02-26 19:42:44
% DurationCPUTime: 0.27s
% Computational Cost: add. (349->58), mult. (870->112), div. (0->0), fcn. (970->14), ass. (0->56)
t325 = sin(pkin(13));
t328 = sin(pkin(6));
t334 = sin(qJ(3));
t336 = cos(qJ(3));
t329 = cos(pkin(13));
t331 = cos(pkin(7));
t356 = t329 * t331;
t327 = sin(pkin(7));
t332 = cos(pkin(6));
t359 = t327 * t332;
t307 = (t325 * t336 + t334 * t356) * t328 + t334 * t359;
t330 = cos(pkin(12));
t326 = sin(pkin(12));
t361 = t326 * t332;
t316 = -t325 * t361 + t329 * t330;
t315 = -t325 * t330 - t329 * t361;
t360 = t327 * t328;
t341 = t315 * t331 + t326 * t360;
t303 = t316 * t336 + t341 * t334;
t355 = t330 * t332;
t314 = t325 * t355 + t326 * t329;
t313 = -t325 * t326 + t329 * t355;
t349 = t330 * t360;
t342 = -t313 * t331 + t349;
t368 = -t314 * t336 + t342 * t334;
t367 = -r_i_i_C(3) - pkin(10) - pkin(9);
t366 = t314 * t334;
t324 = qJ(4) + qJ(5);
t321 = sin(t324);
t323 = qJD(4) + qJD(5);
t363 = t321 * t323;
t322 = cos(t324);
t362 = t322 * t323;
t358 = t328 * t329;
t357 = t328 * t331;
t350 = qJD(3) * t336;
t346 = t331 * t350;
t296 = -t313 * t346 + (t336 * t349 + t366) * qJD(3);
t308 = -t313 * t327 - t330 * t357;
t345 = -t308 * t323 + t296;
t354 = (t345 * t321 + t362 * t368) * r_i_i_C(1) + (t345 * t322 - t363 * t368) * r_i_i_C(2);
t347 = t327 * t350;
t351 = qJD(3) * t334;
t298 = -t326 * t328 * t347 - t315 * t346 + t316 * t351;
t309 = -t315 * t327 + t326 * t357;
t344 = -t309 * t323 + t298;
t353 = (-t303 * t362 + t344 * t321) * r_i_i_C(1) + (t303 * t363 + t344 * t322) * r_i_i_C(2);
t304 = t325 * t328 * t351 - t332 * t347 - t346 * t358;
t312 = -t327 * t358 + t331 * t332;
t343 = -t312 * t323 + t304;
t352 = (-t307 * t362 + t343 * t321) * r_i_i_C(1) + (t307 * t363 + t343 * t322) * r_i_i_C(2);
t335 = cos(qJ(4));
t339 = qJD(3) * (pkin(4) * t335 + r_i_i_C(1) * t322 - r_i_i_C(2) * t321 + pkin(3));
t333 = sin(qJ(4));
t338 = -pkin(4) * qJD(4) * t333 + (-r_i_i_C(1) * t321 - r_i_i_C(2) * t322) * t323;
t1 = [0, 0, t367 * t298 + t338 * (-t316 * t334 + t341 * t336) - t303 * t339 (t298 * t333 + (-t303 * t335 - t309 * t333) * qJD(4)) * pkin(4) + t353, t353, 0; 0, 0, t367 * t296 + t338 * (-t342 * t336 - t366) + t368 * t339 (t296 * t333 + (-t308 * t333 + t335 * t368) * qJD(4)) * pkin(4) + t354, t354, 0; 0, 0, t367 * t304 + t338 * (t336 * t359 + (-t325 * t334 + t336 * t356) * t328) - t307 * t339 (t304 * t333 + (-t307 * t335 - t312 * t333) * qJD(4)) * pkin(4) + t352, t352, 0;];
JaD_transl  = t1;
