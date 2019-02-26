% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:47
% EndTime: 2019-02-26 19:48:47
% DurationCPUTime: 0.37s
% Computational Cost: add. (444->66), mult. (905->120), div. (0->0), fcn. (907->11), ass. (0->48)
t312 = pkin(11) + qJ(4);
t310 = sin(t312);
t311 = cos(t312);
t316 = sin(qJ(6));
t318 = cos(qJ(6));
t333 = r_i_i_C(1) * t318 - r_i_i_C(2) * t316;
t322 = t333 * qJD(6) + qJD(5);
t332 = -r_i_i_C(1) * t316 - r_i_i_C(2) * t318;
t330 = qJ(5) - t332;
t341 = pkin(4) + pkin(9) + r_i_i_C(3);
t320 = (t341 * t310 - t330 * t311) * qJD(4) - t322 * t310;
t313 = sin(pkin(10));
t317 = sin(qJ(2));
t319 = cos(qJ(2));
t347 = cos(pkin(10));
t348 = cos(pkin(6));
t331 = t348 * t347;
t302 = t313 * t319 + t317 * t331;
t314 = sin(pkin(6));
t337 = t314 * t347;
t350 = t302 * t311 - t310 * t337;
t338 = t313 * t348;
t304 = -t317 * t338 + t347 * t319;
t345 = t313 * t314;
t344 = t314 * t317;
t343 = t314 * t319;
t342 = qJD(2) * t317;
t340 = qJD(2) * t343;
t339 = t314 * t342;
t329 = t319 * t331;
t328 = -t304 * t310 + t311 * t345;
t327 = t304 * t311 + t310 * t345;
t326 = pkin(5) + pkin(8) + qJ(3) + t333;
t295 = t310 * t344 - t348 * t311;
t325 = t348 * t310 + t311 * t344;
t324 = -t302 * t310 - t311 * t337;
t323 = t332 * qJD(6) + qJD(3);
t303 = t347 * t317 + t319 * t338;
t321 = -t330 * t310 - t341 * t311 - cos(pkin(11)) * pkin(3) - pkin(2);
t301 = t313 * t317 - t329;
t300 = t304 * qJD(2);
t299 = t303 * qJD(2);
t298 = t302 * qJD(2);
t297 = -qJD(2) * t329 + t313 * t342;
t289 = t325 * qJD(4) + t310 * t340;
t287 = t327 * qJD(4) - t299 * t310;
t285 = t350 * qJD(4) - t297 * t310;
t1 = [0, -t326 * t299 + t321 * t300 + t320 * t303 + t323 * t304, t300, t322 * t327 + t330 * (t328 * qJD(4) - t299 * t311) - t341 * t287, t287 (t287 * t318 - t300 * t316) * r_i_i_C(1) + (-t287 * t316 - t300 * t318) * r_i_i_C(2) + ((-t303 * t318 + t316 * t328) * r_i_i_C(1) + (t303 * t316 + t318 * t328) * r_i_i_C(2)) * qJD(6); 0, -t326 * t297 + t321 * t298 + t320 * t301 + t323 * t302, t298, t322 * t350 + t330 * (t324 * qJD(4) - t297 * t311) - t341 * t285, t285 (t285 * t318 - t298 * t316) * r_i_i_C(1) + (-t285 * t316 - t298 * t318) * r_i_i_C(2) + ((-t301 * t318 + t316 * t324) * r_i_i_C(1) + (t301 * t316 + t318 * t324) * r_i_i_C(2)) * qJD(6); 0 ((t321 * qJD(2) + t323) * t317 + (t326 * qJD(2) - t320) * t319) * t314, t339, t322 * t325 + t330 * (-t295 * qJD(4) + t311 * t340) - t341 * t289, t289 (t289 * t318 - t316 * t339) * r_i_i_C(1) + (-t289 * t316 - t318 * t339) * r_i_i_C(2) + ((-t295 * t316 + t318 * t343) * r_i_i_C(1) + (-t295 * t318 - t316 * t343) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
