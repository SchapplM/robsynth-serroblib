% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP6_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP6_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:51
% EndTime: 2019-02-26 21:48:51
% DurationCPUTime: 0.31s
% Computational Cost: add. (297->77), mult. (888->135), div. (0->0), fcn. (922->10), ass. (0->59)
t316 = cos(pkin(6));
t313 = sin(pkin(11));
t315 = cos(pkin(11));
t321 = cos(qJ(2));
t340 = qJD(2) * t321;
t318 = sin(qJ(2));
t341 = qJD(2) * t318;
t352 = t313 * t341 - t315 * t340;
t294 = t352 * t316;
t328 = t313 * t321 + t315 * t318;
t300 = t328 * t316;
t304 = t313 * t318 - t315 * t321;
t322 = cos(qJ(1));
t319 = sin(qJ(1));
t342 = qJD(1) * t319;
t302 = -t313 * t340 - t315 * t341;
t345 = t319 * t302;
t286 = t345 - t300 * t342 + (-qJD(1) * t304 - t294) * t322;
t314 = sin(pkin(6));
t348 = t314 * t322;
t333 = qJD(4) * t348;
t353 = t286 - t333;
t351 = r_i_i_C(3) + pkin(9);
t350 = pkin(2) * t316;
t349 = t314 * t319;
t347 = t318 * t319;
t346 = t318 * t322;
t344 = t319 * t321;
t343 = t321 * t322;
t288 = t300 * t322 - t304 * t319;
t339 = qJD(4) * t288;
t338 = pkin(2) * t341;
t337 = t314 * t342;
t336 = qJD(1) * t348;
t317 = sin(qJ(4));
t320 = cos(qJ(4));
t332 = r_i_i_C(1) * t317 + r_i_i_C(2) * t320;
t299 = t304 * t316;
t331 = -t299 * t322 - t319 * t328;
t330 = t299 * t319 - t322 * t328;
t329 = t300 * t319 + t304 * t322;
t327 = r_i_i_C(1) * t320 - r_i_i_C(2) * t317 + pkin(3);
t326 = qJD(4) * t332;
t325 = t337 - t339;
t324 = qJD(2) * t328;
t323 = t304 * qJD(2);
t283 = -qJD(1) * t288 + t319 * t294 + t322 * t302;
t312 = pkin(2) * t321 + pkin(1);
t308 = t320 * t333;
t303 = -qJD(3) * t314 + t340 * t350;
t301 = t318 * t350 + (-pkin(8) - qJ(3)) * t314;
t298 = t328 * t314;
t295 = t316 * t324;
t292 = t352 * t314;
t285 = qJD(1) * t330 - t322 * t295 + t319 * t323;
t282 = qJD(1) * t331 - t295 * t319 - t322 * t323;
t280 = t317 * t336 + t283 * t320 + (t317 * t329 + t320 * t349) * qJD(4);
t279 = t320 * t336 - t283 * t317 + (-t317 * t349 + t320 * t329) * qJD(4);
t1 = [(-t286 * t320 + t317 * t339 + t308) * r_i_i_C(1) + (t317 * t353 + t320 * t339) * r_i_i_C(2) - t286 * pkin(3) + t319 * t338 - t322 * t303 + t351 * t285 + (-t322 * t312 + (-t314 * t332 + t301) * t319) * qJD(1), t351 * t283 - t330 * t326 - t327 * t282 + ((t316 * t347 - t343) * qJD(2) + (-t316 * t343 + t347) * qJD(1)) * pkin(2), t336, r_i_i_C(1) * t279 - r_i_i_C(2) * t280, 0, 0; -t322 * t338 + t283 * pkin(3) + t280 * r_i_i_C(1) + t279 * r_i_i_C(2) - t319 * t303 + t351 * t282 + (-t301 * t322 - t312 * t319) * qJD(1), -t351 * (qJD(1) * t329 + t322 * t294 - t345) - t331 * t326 + t327 * t285 + ((-t316 * t346 - t344) * qJD(2) + (-t316 * t344 - t346) * qJD(1)) * pkin(2), t337, t308 * r_i_i_C(2) + (r_i_i_C(1) * t325 - t286 * r_i_i_C(2)) * t320 + (-r_i_i_C(1) * t353 - t325 * r_i_i_C(2)) * t317, 0, 0; 0, -t351 * t292 + (t304 * t326 - t324 * t327 - t338) * t314, 0, t332 * t292 + ((-t298 * t320 - t316 * t317) * r_i_i_C(1) + (t298 * t317 - t316 * t320) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
