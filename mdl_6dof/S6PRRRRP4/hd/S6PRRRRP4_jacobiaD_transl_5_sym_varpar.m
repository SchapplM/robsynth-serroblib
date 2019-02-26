% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRP4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:17:00
% EndTime: 2019-02-26 20:17:00
% DurationCPUTime: 0.42s
% Computational Cost: add. (487->86), mult. (1093->149), div. (0->0), fcn. (1098->12), ass. (0->63)
t334 = sin(qJ(3));
t337 = cos(qJ(3));
t336 = cos(qJ(4));
t330 = qJ(4) + qJ(5);
t327 = sin(t330);
t328 = cos(t330);
t353 = t328 * r_i_i_C(1) - t327 * r_i_i_C(2);
t349 = t336 * pkin(4) + pkin(3) + t353;
t378 = r_i_i_C(3) + pkin(10) + pkin(9);
t384 = (t334 * t349 - t337 * t378) * qJD(3);
t335 = sin(qJ(2));
t338 = cos(qJ(2));
t331 = sin(pkin(11));
t377 = cos(pkin(6));
t362 = t331 * t377;
t376 = cos(pkin(11));
t320 = -t335 * t362 + t376 * t338;
t352 = t327 * r_i_i_C(1) + t328 * r_i_i_C(2);
t329 = qJD(4) + qJD(5);
t333 = sin(qJ(4));
t379 = t333 * pkin(4);
t383 = qJD(4) * t379 + t329 * t352;
t375 = t327 * t329;
t374 = t328 * t329;
t332 = sin(pkin(6));
t373 = t332 * t334;
t372 = t332 * t337;
t371 = t332 * t338;
t350 = t377 * t376;
t318 = t331 * t338 + t335 * t350;
t314 = t318 * qJD(2);
t361 = t332 * t376;
t343 = -t318 * t337 + t334 * t361;
t357 = -t329 * t343 - t314;
t348 = t338 * t350;
t367 = qJD(2) * t335;
t313 = -qJD(2) * t348 + t331 * t367;
t345 = -t318 * t334 - t337 * t361;
t304 = qJD(3) * t345 - t313 * t337;
t317 = t331 * t335 - t348;
t359 = -t317 * t329 - t304;
t370 = (t327 * t359 - t328 * t357) * r_i_i_C(1) + (t327 * t357 + t328 * t359) * r_i_i_C(2);
t310 = t320 * t337 + t331 * t373;
t316 = t320 * qJD(2);
t356 = t310 * t329 - t316;
t319 = t335 * t376 + t338 * t362;
t315 = t319 * qJD(2);
t347 = -t320 * t334 + t331 * t372;
t306 = qJD(3) * t347 - t315 * t337;
t358 = -t319 * t329 - t306;
t369 = (t327 * t358 - t328 * t356) * r_i_i_C(1) + (t327 * t356 + t328 * t358) * r_i_i_C(2);
t322 = t334 * t377 + t335 * t372;
t364 = t332 * t367;
t346 = -t322 * t329 + t364;
t344 = -t335 * t373 + t337 * t377;
t363 = qJD(2) * t371;
t312 = qJD(3) * t344 + t337 * t363;
t351 = t329 * t371 - t312;
t368 = (t327 * t351 + t328 * t346) * r_i_i_C(1) + (-t327 * t346 + t328 * t351) * r_i_i_C(2);
t366 = qJD(4) * t336;
t341 = -t334 * t378 - t337 * t349 - pkin(2);
t340 = t383 * t337 + t384;
t1 = [0 (-t315 * t327 + t320 * t374) * r_i_i_C(1) + (-t315 * t328 - t320 * t375) * r_i_i_C(2) - t315 * pkin(8) + (-t315 * t333 + t320 * t366) * pkin(4) + t341 * t316 + t340 * t319, t378 * t306 - t383 * t347 + t349 * (-qJD(3) * t310 + t315 * t334) (-t306 * t333 + t316 * t336 + (-t310 * t336 - t319 * t333) * qJD(4)) * pkin(4) + t369, t369, 0; 0 (-t313 * t327 + t318 * t374) * r_i_i_C(1) + (-t313 * t328 - t318 * t375) * r_i_i_C(2) - t313 * pkin(8) + (-t313 * t333 + t318 * t366) * pkin(4) + t341 * t314 + t340 * t317, t378 * t304 - t383 * t345 + t349 * (qJD(3) * t343 + t313 * t334) (-t304 * t333 + t314 * t336 + (-t317 * t333 + t336 * t343) * qJD(4)) * pkin(4) + t370, t370, 0; 0 ((pkin(4) * t366 + qJD(2) * t341 + t329 * t353) * t335 + (qJD(2) * pkin(8) + (-qJD(4) * t337 + qJD(2)) * t379 - t384 + t352 * (-t329 * t337 + qJD(2))) * t338) * t332, t378 * t312 - t383 * t344 + t349 * (-qJD(3) * t322 - t334 * t363) (t336 * t364 - t312 * t333 + (-t322 * t336 + t333 * t371) * qJD(4)) * pkin(4) + t368, t368, 0;];
JaD_transl  = t1;
