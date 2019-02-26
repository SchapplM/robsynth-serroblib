% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR8_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR8_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobiaD_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:18
% EndTime: 2019-02-26 20:08:18
% DurationCPUTime: 0.22s
% Computational Cost: add. (224->62), mult. (748->112), div. (0->0), fcn. (768->10), ass. (0->46)
t309 = sin(pkin(7));
t330 = (pkin(9) + r_i_i_C(1)) * t309;
t308 = sin(pkin(12));
t311 = cos(pkin(12));
t315 = sin(qJ(2));
t313 = cos(pkin(6));
t317 = cos(qJ(2));
t335 = t313 * t317;
t302 = -t308 * t315 + t311 * t335;
t314 = sin(qJ(3));
t316 = cos(qJ(3));
t336 = t313 * t315;
t322 = t308 * t336 - t311 * t317;
t312 = cos(pkin(7));
t323 = t308 * t335 + t311 * t315;
t310 = sin(pkin(6));
t340 = t309 * t310;
t324 = t308 * t340 - t312 * t323;
t348 = t324 * t314 - t316 * t322;
t303 = t308 * t317 + t311 * t336;
t326 = -t302 * t312 + t311 * t340;
t347 = -t303 * t316 + t326 * t314;
t346 = pkin(3) - r_i_i_C(2);
t344 = r_i_i_C(3) + qJ(4);
t339 = t309 * t313;
t338 = t312 * t314;
t337 = t312 * t316;
t334 = t314 * t315;
t333 = t314 * t317;
t332 = t315 * t316;
t331 = t316 * t317;
t328 = qJD(3) * t339;
t327 = -t302 * t314 - t303 * t337;
t325 = t314 * t323 + t322 * t337;
t321 = t312 * t331 - t334;
t320 = t312 * t332 + t333;
t319 = t312 * t333 + t332;
t318 = t312 * t334 - t331;
t301 = t322 * qJD(2);
t300 = t323 * qJD(2);
t299 = t303 * qJD(2);
t298 = t302 * qJD(2);
t294 = t314 * t328 + (t320 * qJD(2) + t319 * qJD(3)) * t310;
t288 = t348 * qJD(3) - t300 * t314 - t301 * t337;
t286 = -t347 * qJD(3) + t298 * t314 + t299 * t337;
t1 = [0, -t325 * qJD(4) + t301 * pkin(2) - t300 * t330 + t346 * (t325 * qJD(3) + t300 * t338 + t301 * t316) + t344 * (-t300 * t337 + t301 * t314 + (-t316 * t323 + t322 * t338) * qJD(3)) t348 * qJD(4) + t344 * (t301 * t338 - t300 * t316 + (t314 * t322 + t324 * t316) * qJD(3)) - t346 * t288, t288, 0, 0; 0, -t327 * qJD(4) - t299 * pkin(2) + t298 * t330 + t346 * (t327 * qJD(3) - t298 * t338 - t299 * t316) + t344 * (t298 * t337 - t299 * t314 + (t302 * t316 - t303 * t338) * qJD(3)) -t347 * qJD(4) + t344 * (-t299 * t338 + t298 * t316 + (-t303 * t314 - t326 * t316) * qJD(3)) - t346 * t286, t286, 0, 0; 0 (-t346 * (t319 * qJD(2) + t320 * qJD(3)) - t344 * (-t321 * qJD(2) + t318 * qJD(3)) + t320 * qJD(4) + (-t315 * pkin(2) + t317 * t330) * qJD(2)) * t310 -(-t319 * t310 - t314 * t339) * qJD(4) + t344 * (t316 * t328 + (-t318 * qJD(2) + t321 * qJD(3)) * t310) - t346 * t294, t294, 0, 0;];
JaD_transl  = t1;
