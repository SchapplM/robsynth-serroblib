% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR13_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR13_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobiaD_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:40
% EndTime: 2019-02-26 20:55:40
% DurationCPUTime: 0.24s
% Computational Cost: add. (232->51), mult. (778->91), div. (0->0), fcn. (800->10), ass. (0->51)
t324 = cos(pkin(12));
t326 = cos(pkin(6));
t321 = sin(pkin(12));
t328 = sin(qJ(1));
t349 = t328 * t321;
t341 = t326 * t349;
t330 = cos(qJ(1));
t345 = qJD(1) * t330;
t313 = -qJD(1) * t341 + t324 * t345;
t327 = sin(qJ(3));
t329 = cos(qJ(3));
t347 = t330 * t321;
t348 = t328 * t324;
t315 = t326 * t347 + t348;
t346 = t330 * t324;
t314 = -t326 * t346 + t349;
t322 = sin(pkin(7));
t325 = cos(pkin(7));
t323 = sin(pkin(6));
t351 = t323 * t330;
t337 = t314 * t325 + t322 * t351;
t331 = -t315 * t329 + t337 * t327;
t334 = t326 * t348 + t347;
t312 = t334 * qJD(1);
t352 = t323 * t328;
t340 = qJD(1) * t352;
t333 = -t312 * t325 + t322 * t340;
t302 = -t331 * qJD(3) + t313 * t327 - t333 * t329;
t332 = t315 * t327 + t337 * t329;
t370 = -t332 * qJD(3) + t313 * t329 + t333 * t327;
t362 = r_i_i_C(1) + pkin(9);
t367 = t362 * t325 + qJ(2);
t353 = t322 * t326;
t366 = (t324 * t325 * t327 + t321 * t329) * t323 + t327 * t353;
t317 = -t341 + t346;
t336 = t322 * t352 - t334 * t325;
t364 = t317 * t329 + t336 * t327;
t361 = -r_i_i_C(2) + pkin(3);
t360 = r_i_i_C(3) + qJ(4);
t355 = t317 * t327;
t350 = t325 * t329;
t344 = t323 * qJD(2);
t343 = t362 * t322;
t339 = t323 * t345;
t338 = t322 * t339;
t311 = t315 * qJD(1);
t310 = t314 * qJD(1);
t307 = t366 * qJD(3);
t301 = -qJD(3) * t355 + (t310 * t325 + t338) * t327 + (t336 * qJD(3) - t311) * t329;
t300 = t364 * qJD(3) - t310 * t350 - t311 * t327 - t329 * t338;
t1 = [-t332 * qJD(4) - t313 * pkin(2) + t330 * t344 - t312 * t343 - t361 * t370 - t360 * t302 + (-t330 * pkin(1) - t367 * t352) * qJD(1), t339, t364 * qJD(4) - t361 * t300 + t360 * t301, t300, 0, 0; -(t336 * t329 - t355) * qJD(4) - t311 * pkin(2) + t328 * t344 - t310 * t343 + t361 * t301 + t360 * t300 + (-t328 * pkin(1) + t367 * t351) * qJD(1), t340, -t331 * qJD(4) - t361 * t302 + t360 * t370, t302, 0, 0; 0, 0, t366 * qJD(4) - t361 * t307 - t360 * (-t329 * t353 + (t321 * t327 - t324 * t350) * t323) * qJD(3), t307, 0, 0;];
JaD_transl  = t1;
