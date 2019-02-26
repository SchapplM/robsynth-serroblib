% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:47
% EndTime: 2019-02-26 22:30:47
% DurationCPUTime: 0.40s
% Computational Cost: add. (792->76), mult. (564->97), div. (0->0), fcn. (413->12), ass. (0->65)
t283 = qJ(2) + qJ(3);
t281 = qJ(4) + t283;
t274 = pkin(11) + t281;
t272 = sin(t274);
t287 = cos(qJ(6));
t339 = r_i_i_C(1) * t287 + pkin(5);
t304 = t339 * t272;
t284 = sin(qJ(6));
t322 = qJD(6) * t272;
t273 = cos(t274);
t282 = qJD(2) + qJD(3);
t277 = qJD(4) + t282;
t327 = t273 * t277;
t342 = t284 * t327 + t287 * t322;
t336 = pkin(10) + r_i_i_C(3);
t316 = t336 * t273;
t319 = r_i_i_C(2) * t272 * t284;
t341 = t316 + t319;
t275 = sin(t281);
t285 = sin(qJ(2));
t329 = pkin(2) * qJD(2);
t318 = t285 * t329;
t278 = sin(t283);
t335 = pkin(3) * t282;
t320 = t278 * t335;
t332 = pkin(4) * t277;
t340 = -t275 * t332 + (-pkin(5) * t272 + t316) * t277 - t318 - t320;
t311 = t284 * t322;
t337 = r_i_i_C(1) * t311 + t342 * r_i_i_C(2);
t276 = cos(t281);
t279 = cos(t283);
t317 = t336 * t272;
t297 = -t273 * t339 - t317;
t293 = -t276 * t332 + t297 * t277 - t279 * t335;
t334 = pkin(4) * t275;
t333 = pkin(4) * t276;
t328 = t272 * t277;
t326 = t277 * t287;
t289 = cos(qJ(1));
t325 = t287 * t289;
t286 = sin(qJ(1));
t324 = qJD(1) * t286;
t323 = qJD(1) * t289;
t321 = qJD(6) * t273;
t306 = -qJD(1) + t321;
t305 = qJD(1) * t273 - qJD(6);
t266 = -pkin(3) * t278 - t334;
t303 = t337 * t289 + t324 * t304;
t302 = t306 * t284;
t301 = t337 * t286 + t341 * t323;
t288 = cos(qJ(2));
t299 = -pkin(2) * t288 - pkin(3) * t279 - pkin(5) * t273 - pkin(1) - t317 - t333;
t298 = -t304 - t334;
t296 = t305 * t286 + t289 * t328;
t294 = -t288 * t329 + t293;
t292 = t277 * (t297 - t333);
t291 = (-r_i_i_C(1) * t284 - r_i_i_C(2) * t287) * t321 + t336 * t327 + (t319 + t298) * t277;
t290 = t291 - t320;
t280 = -qJ(5) - pkin(9) - pkin(8) - pkin(7);
t258 = -pkin(2) * t285 + t266;
t251 = -t305 * t325 + (t272 * t326 + t302) * t286;
t250 = t306 * t287 * t286 + (-t286 * t328 + t305 * t289) * t284;
t249 = t296 * t287 + t289 * t302;
t248 = t296 * t284 - t306 * t325;
t1 = [t251 * r_i_i_C(1) + t250 * r_i_i_C(2) + t289 * qJD(5) - t340 * t286 + (t280 * t286 + t299 * t289) * qJD(1) (-t258 - t341) * t324 + t294 * t289 + t303 (-t266 - t341) * t324 + t293 * t289 + t303 (-t341 + t334) * t324 + t289 * t292 + t303, t323, r_i_i_C(1) * t248 + r_i_i_C(2) * t249; -t249 * r_i_i_C(1) + t248 * r_i_i_C(2) + t286 * qJD(5) + t340 * t289 + (-t280 * t289 + t299 * t286) * qJD(1) (t258 - t304) * t323 + t294 * t286 + t301 (t266 - t304) * t323 + t293 * t286 + t301, t286 * t292 + t298 * t323 + t301, t324, -r_i_i_C(1) * t250 + r_i_i_C(2) * t251; 0, t290 - t318, t290, t291, 0 (-t273 * t326 + t311) * r_i_i_C(2) - t342 * r_i_i_C(1);];
JaD_transl  = t1;
