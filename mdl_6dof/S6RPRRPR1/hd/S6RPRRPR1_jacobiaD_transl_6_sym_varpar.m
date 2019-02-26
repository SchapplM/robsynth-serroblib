% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:00:51
% EndTime: 2019-02-26 21:00:51
% DurationCPUTime: 0.31s
% Computational Cost: add. (532->65), mult. (430->88), div. (0->0), fcn. (321->12), ass. (0->57)
t281 = qJ(3) + qJ(4);
t275 = pkin(11) + t281;
t271 = sin(t275);
t284 = cos(qJ(6));
t333 = r_i_i_C(1) * t284 + pkin(5);
t297 = t333 * t271;
t272 = cos(t275);
t315 = qJD(6) * t271;
t279 = qJD(3) + qJD(4);
t282 = sin(qJ(6));
t319 = t279 * t282;
t336 = t272 * t319 + t284 * t315;
t327 = pkin(9) + r_i_i_C(3);
t335 = t327 * t272;
t293 = -r_i_i_C(2) * t271 * t282 - t335;
t276 = sin(t281);
t283 = sin(qJ(3));
t321 = pkin(3) * qJD(3);
t312 = t283 * t321;
t324 = pkin(4) * t279;
t334 = -t276 * t324 + (-pkin(5) * t271 + t335) * t279 - t312;
t304 = t282 * t315;
t330 = r_i_i_C(1) * t304 + t336 * r_i_i_C(2);
t298 = qJD(1) * t272 - qJD(6);
t329 = t284 * t298;
t314 = qJD(6) * t272;
t299 = -qJD(1) + t314;
t309 = t271 * t319;
t328 = t299 * t284 - t309;
t326 = pkin(4) * t276;
t277 = cos(t281);
t325 = pkin(4) * t277;
t318 = t279 * t284;
t280 = qJ(1) + pkin(10);
t273 = sin(t280);
t317 = qJD(1) * t273;
t274 = cos(t280);
t316 = qJD(1) * t274;
t311 = t327 * t271;
t296 = t330 * t274 + t317 * t297;
t295 = t298 * t282;
t294 = t330 * t273 - t293 * t316;
t285 = cos(qJ(3));
t292 = -pkin(3) * t285 - pkin(5) * t272 - pkin(2) - t311 - t325;
t291 = -t297 - t326;
t290 = -t272 * t333 - t311;
t289 = t271 * t318 + t299 * t282;
t288 = -t277 * t324 + t290 * t279 - t285 * t321;
t287 = (t290 - t325) * t279;
t286 = r_i_i_C(2) * t309 + (-r_i_i_C(1) * t282 - r_i_i_C(2) * t284) * t314 + (t291 + t335) * t279;
t278 = -qJ(5) - pkin(8) - pkin(7);
t270 = -pkin(3) * t283 - t326;
t252 = t289 * t273 - t274 * t329;
t251 = t328 * t273 + t274 * t295;
t250 = t273 * t329 + t289 * t274;
t249 = t273 * t295 - t328 * t274;
t1 = [t252 * r_i_i_C(1) + t251 * r_i_i_C(2) + t274 * qJD(5) - t334 * t273 + (-cos(qJ(1)) * pkin(1) + t273 * t278 + t292 * t274) * qJD(1), 0 (-t270 + t293) * t317 + t288 * t274 + t296 (t293 + t326) * t317 + t274 * t287 + t296, t316, r_i_i_C(1) * t249 + r_i_i_C(2) * t250; -t250 * r_i_i_C(1) + t249 * r_i_i_C(2) + t273 * qJD(5) + t334 * t274 + (-sin(qJ(1)) * pkin(1) - t274 * t278 + t292 * t273) * qJD(1), 0 (t270 - t297) * t316 + t288 * t273 + t294, t273 * t287 + t291 * t316 + t294, t317, -r_i_i_C(1) * t251 + r_i_i_C(2) * t252; 0, 0, t286 - t312, t286, 0 (-t272 * t318 + t304) * r_i_i_C(2) - t336 * r_i_i_C(1);];
JaD_transl  = t1;
