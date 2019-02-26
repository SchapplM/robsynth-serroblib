% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRP5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:52:32
% EndTime: 2019-02-26 19:52:33
% DurationCPUTime: 0.41s
% Computational Cost: add. (331->78), mult. (1015->134), div. (0->0), fcn. (1016->10), ass. (0->61)
t312 = sin(qJ(5));
t315 = cos(qJ(5));
t351 = pkin(5) + r_i_i_C(1);
t354 = -t312 * r_i_i_C(2) + t351 * t315;
t322 = t354 * qJD(5);
t313 = sin(qJ(4));
t316 = cos(qJ(4));
t328 = pkin(4) + t354;
t349 = r_i_i_C(3) + qJ(6) + pkin(9);
t352 = t328 * t313 - t349 * t316 + qJ(3);
t348 = cos(pkin(6));
t308 = sin(pkin(10));
t310 = cos(pkin(10));
t314 = sin(qJ(2));
t317 = cos(qJ(2));
t336 = t317 * t348;
t297 = t308 * t336 + t310 * t314;
t347 = t297 * t316;
t309 = sin(pkin(6));
t346 = t309 * t313;
t345 = t309 * t314;
t344 = t309 * t316;
t343 = t309 * t317;
t342 = qJD(2) * t314;
t341 = qJD(2) * t317;
t340 = t309 * t341;
t339 = qJD(4) * t346;
t338 = t309 * t342;
t337 = t314 * t348;
t335 = t348 * t316;
t334 = t308 * t337;
t333 = t310 * t336;
t294 = -qJD(2) * t334 + t310 * t341;
t277 = -qJD(4) * t347 - t294 * t313 + t308 * t339;
t293 = t297 * qJD(2);
t332 = t277 * t312 - t293 * t315;
t296 = t308 * t317 + t310 * t337;
t292 = t296 * qJD(2);
t295 = t308 * t314 - t333;
t326 = t295 * t316 + t310 * t346;
t279 = t326 * qJD(4) + t292 * t313;
t291 = -qJD(2) * t333 + t308 * t342;
t331 = -t279 * t312 - t291 * t315;
t282 = t297 * t313 + t308 * t344;
t298 = t310 * t317 - t334;
t330 = -t282 * t315 - t298 * t312;
t284 = -t295 * t313 + t310 * t344;
t329 = t284 * t315 - t296 * t312;
t327 = -t315 * r_i_i_C(2) - t351 * t312;
t300 = -t313 * t343 + t335;
t325 = -t300 * t315 - t312 * t345;
t324 = -pkin(2) - pkin(8) + t327;
t320 = t348 * t313 + t316 * t343;
t285 = t320 * qJD(4) - t313 * t338;
t323 = t285 * t312 + t315 * t340;
t321 = qJD(5) * t327;
t318 = -t316 * qJD(6) + qJD(3) + t313 * t321 + (t349 * t313 + t328 * t316) * qJD(4);
t286 = qJD(4) * t335 - t316 * t338 - t317 * t339;
t280 = t284 * qJD(4) + t292 * t316;
t278 = t282 * qJD(4) - t294 * t316;
t1 = [0, -t293 * t352 + t324 * t294 - t297 * t322 + t318 * t298, t294, t282 * qJD(6) - t349 * t277 - t328 * t278 + (-t308 * t346 + t347) * t321, t332 * r_i_i_C(1) + (t277 * t315 + t293 * t312) * r_i_i_C(2) + (t330 * r_i_i_C(1) + (t282 * t312 - t298 * t315) * r_i_i_C(2)) * qJD(5) + (t330 * qJD(5) + t332) * pkin(5), t278; 0, -t291 * t352 + t324 * t292 - t295 * t322 + t318 * t296, t292, -t284 * qJD(6) + t349 * t279 + t328 * t280 + t326 * t321, t331 * r_i_i_C(1) + (-t279 * t315 + t291 * t312) * r_i_i_C(2) + (t329 * r_i_i_C(1) + (-t284 * t312 - t296 * t315) * r_i_i_C(2)) * qJD(5) + (t329 * qJD(5) + t331) * pkin(5), -t280; 0 ((t352 * qJD(2) + t322) * t317 + (t324 * qJD(2) + t318) * t314) * t309, t338, t300 * qJD(6) - t349 * t285 - t328 * t286 - t320 * t321, t323 * r_i_i_C(1) + (t285 * t315 - t312 * t340) * r_i_i_C(2) + (t325 * r_i_i_C(1) + (t300 * t312 - t315 * t345) * r_i_i_C(2)) * qJD(5) + (t325 * qJD(5) + t323) * pkin(5), t286;];
JaD_transl  = t1;
