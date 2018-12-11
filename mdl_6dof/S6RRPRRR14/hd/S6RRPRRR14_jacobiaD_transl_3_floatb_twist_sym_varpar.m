% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRR14
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR14_jacobiaD_transl_3_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_3_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_3_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobiaD_transl_3_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_transl_3_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:22
% EndTime: 2018-12-10 18:38:23
% DurationCPUTime: 0.22s
% Computational Cost: add. (407->69), mult. (584->110), div. (0->0), fcn. (432->17), ass. (0->48)
t357 = pkin(6) + qJ(2);
t331 = cos(t357) / 0.2e1;
t358 = pkin(6) - qJ(2);
t335 = cos(t358);
t370 = t331 - t335 / 0.2e1;
t323 = t335 / 0.2e1 + t331;
t343 = sin(qJ(2));
t344 = sin(qJ(1));
t346 = cos(qJ(1));
t369 = -t346 * t323 + t344 * t343;
t316 = t323 * qJD(2);
t354 = -sin(t358) / 0.2e1;
t355 = sin(t357);
t322 = t355 / 0.2e1 + t354;
t345 = cos(qJ(2));
t353 = t346 * t322 + t344 * t345;
t359 = qJD(2) * t346;
t309 = t353 * qJD(1) + t344 * t316 + t343 * t359;
t366 = -r_i_i_C(3) - qJ(3);
t348 = t322 * qJD(2);
t308 = t369 * qJD(1) + t344 * t348 - t345 * t359;
t339 = sin(pkin(7));
t365 = t308 * t339;
t340 = sin(pkin(6));
t364 = t340 * t344;
t363 = t340 * t346;
t360 = qJD(2) * t344;
t342 = cos(pkin(7));
t356 = qJD(1) * t340 * t342;
t352 = t344 * t322 - t346 * t345;
t351 = t344 * t323 + t346 * t343;
t336 = pkin(7) + pkin(14);
t327 = sin(t336) / 0.2e1;
t337 = pkin(7) - pkin(14);
t328 = cos(t337) / 0.2e1;
t332 = sin(t337);
t333 = cos(t336);
t350 = -r_i_i_C(1) * (t328 - t333 / 0.2e1) - r_i_i_C(2) * (t327 + t332 / 0.2e1) - pkin(10);
t347 = t352 * qJD(1) - t346 * t316 + t343 * t360;
t341 = cos(pkin(14));
t338 = sin(pkin(14));
t320 = t328 + t333 / 0.2e1;
t319 = t327 - t332 / 0.2e1;
t317 = t370 * qJD(2);
t315 = (-t355 / 0.2e1 + t354) * qJD(2);
t311 = t351 * qJD(1) + t345 * t360 + t346 * t348;
t307 = t346 * t356 - t365;
t1 = [(t311 * t319 + t341 * t347) * r_i_i_C(1) + (t311 * t320 - t338 * t347) * r_i_i_C(2) + t347 * pkin(2) + t342 * qJD(3) * t363 + (-t369 * qJD(3) + t366 * t311) * t339 + (-t346 * pkin(1) + (t366 * t342 + t350) * t364) * qJD(1) (t308 * t341 + t309 * t319) * r_i_i_C(1) + (-t308 * t338 + t309 * t320) * r_i_i_C(2) + t308 * pkin(2) + (-t352 * qJD(3) + t366 * t309) * t339, t307, 0, 0, 0; (t308 * t319 - t309 * t341) * r_i_i_C(1) + (t308 * t320 + t309 * t338) * r_i_i_C(2) + t307 * r_i_i_C(3) - t309 * pkin(2) - qJ(3) * t365 + (t351 * t339 + t342 * t364) * qJD(3) + (-t344 * pkin(1) + (qJ(3) * t342 - t350) * t363) * qJD(1) (-t311 * t341 + t319 * t347) * r_i_i_C(1) + (t311 * t338 + t320 * t347) * r_i_i_C(2) - t311 * pkin(2) + (t353 * qJD(3) + t347 * t366) * t339, t311 * t339 + t344 * t356, 0, 0, 0; 0 (t315 * t319 + t317 * t341) * r_i_i_C(1) + (t315 * t320 - t317 * t338) * r_i_i_C(2) + t317 * pkin(2) + (-t370 * qJD(3) + t366 * t315) * t339, -t317 * t339, 0, 0, 0;];
JaD_transl  = t1;
