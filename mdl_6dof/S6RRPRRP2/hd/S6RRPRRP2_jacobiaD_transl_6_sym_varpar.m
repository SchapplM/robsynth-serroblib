% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:46:36
% EndTime: 2019-02-26 21:46:37
% DurationCPUTime: 0.36s
% Computational Cost: add. (677->75), mult. (720->100), div. (0->0), fcn. (578->10), ass. (0->56)
t326 = sin(qJ(5));
t329 = cos(qJ(5));
t331 = cos(qJ(1));
t367 = qJD(5) * t331;
t328 = sin(qJ(1));
t370 = qJD(1) * t328;
t339 = t326 * t367 + t329 * t370;
t377 = r_i_i_C(1) + pkin(5);
t376 = r_i_i_C(3) + qJ(6);
t382 = t376 * t326;
t325 = qJ(2) + pkin(10);
t322 = qJ(4) + t325;
t319 = cos(t322);
t378 = pkin(9) + r_i_i_C(2);
t362 = t378 * t319;
t308 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t325);
t298 = t308 * qJD(2);
t324 = qJD(2) + qJD(4);
t361 = t378 * t324;
t366 = qJD(6) * t326;
t318 = sin(t322);
t375 = t318 * t324;
t381 = -pkin(4) * t375 + (t361 + t366) * t319 + t298;
t368 = qJD(5) * t329;
t337 = -t376 * t368 - t366;
t374 = t319 * t324;
t373 = t319 * t328;
t372 = t328 * t326;
t371 = t331 * t329;
t369 = qJD(1) * t331;
t365 = t329 * qJD(6);
t364 = t377 * t326;
t363 = t378 * t318;
t360 = t328 * t375;
t359 = t331 * t375;
t354 = qJD(5) * t372;
t352 = t329 * t367;
t351 = qJD(3) - t365;
t350 = t377 * t318 * t354 + t369 * t362;
t344 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t325);
t343 = t319 * t371 + t372;
t342 = -pkin(4) * t319 - pkin(1) + t344 - t363;
t341 = (t377 * t339 + (pkin(4) + t382) * t370) * t318;
t340 = -t377 * t329 - t382;
t338 = t326 * t369 + t328 * t368;
t336 = -pkin(4) + t340;
t335 = t336 * t324;
t334 = t318 * t335 + t378 * t374 + (-qJD(5) * t364 - t337) * t319;
t333 = t337 * t318 + (t336 * t319 - t363) * t324;
t332 = t344 * qJD(2) + t333;
t323 = -pkin(8) - qJ(3) - pkin(7);
t287 = t343 * qJD(1) - t319 * t354 - t329 * t360 - t352;
t286 = t338 * t319 - t326 * t360 - t339;
t285 = t339 * t319 + t329 * t359 - t338;
t284 = t326 * t359 - t319 * t352 - t354 + (t319 * t372 + t371) * qJD(1);
t1 = [t351 * t331 - t377 * t287 - t376 * t286 - t381 * t328 + (t328 * t323 + t342 * t331) * qJD(1) (-t308 - t362) * t370 + t332 * t331 + t341, t369, t333 * t331 - t362 * t370 + t341, t343 * qJD(6) + t377 * t284 - t376 * t285, -t284; t351 * t328 - t377 * t285 - t376 * t284 + t381 * t331 + (-t331 * t323 + t342 * t328) * qJD(1) (t336 * t318 + t308) * t369 + t332 * t328 + t350, t370, t335 * t373 + ((-t361 + t337) * t328 + t336 * t369) * t318 + t350 -(t331 * t326 - t329 * t373) * qJD(6) + t376 * t287 - t377 * t286, t286; 0, t298 + t334, 0, t334 (t376 * t329 - t364) * t374 + (t340 * qJD(5) + t365) * t318, t318 * t368 + t326 * t374;];
JaD_transl  = t1;
