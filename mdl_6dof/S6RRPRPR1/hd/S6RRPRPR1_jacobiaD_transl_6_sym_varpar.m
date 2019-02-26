% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:51
% EndTime: 2019-02-26 21:37:51
% DurationCPUTime: 0.31s
% Computational Cost: add. (539->73), mult. (455->102), div. (0->0), fcn. (351->12), ass. (0->55)
t276 = pkin(11) + qJ(6);
t270 = sin(t276);
t272 = cos(t276);
t278 = qJ(2) + pkin(10);
t274 = qJ(4) + t278;
t267 = sin(t274);
t312 = qJD(6) * t267;
t268 = cos(t274);
t277 = qJD(2) + qJD(4);
t319 = t268 * t277;
t332 = t270 * t319 + t272 * t312;
t269 = cos(pkin(11)) * pkin(5) + pkin(4);
t331 = r_i_i_C(1) * t272 + t269;
t262 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t278);
t257 = t262 * qJD(2);
t266 = t267 * qJD(5);
t280 = -pkin(9) - qJ(5);
t321 = r_i_i_C(3) - t280;
t330 = (-t267 * t269 + t321 * t268) * t277 + (pkin(5) * sin(pkin(11)) + pkin(7) + pkin(8) + qJ(3)) * qJD(1) + t257 + t266;
t305 = t270 * t312;
t329 = r_i_i_C(1) * t305 + t332 * r_i_i_C(2) + qJD(5) * t268;
t284 = cos(qJ(1));
t297 = qJD(6) * t268 - qJD(1);
t328 = t284 * t297;
t296 = qJD(1) * t268 - qJD(6);
t282 = sin(qJ(1));
t317 = t277 * t282;
t310 = t267 * t317;
t326 = t296 * t284 - t310;
t323 = r_i_i_C(2) * t270;
t322 = r_i_i_C(3) * t268;
t318 = t268 * t280;
t316 = t277 * t284;
t315 = qJD(1) * t282;
t314 = qJD(1) * t284;
t311 = t267 * t323;
t309 = t267 * t316;
t307 = t267 * t315;
t306 = t267 * t314;
t295 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t278);
t294 = t331 * t277;
t293 = t297 * t282;
t292 = t280 * t310 + t329 * t282 + t306 * t323 + t314 * t322;
t291 = -t267 * t331 - t318;
t290 = t280 * t309 + t329 * t284 + t331 * t307 + t315 * t318;
t289 = (-r_i_i_C(3) * t267 - t268 * t331) * t277;
t288 = t296 * t282 + t309;
t287 = qJD(3) + (-t321 * t267 - t268 * t269 - pkin(1) + t295) * qJD(1);
t286 = t295 * qJD(2) + t289;
t285 = t277 * t311 + r_i_i_C(3) * t319 + t266 - t267 * t294 + (-t277 * t280 + (-r_i_i_C(1) * t270 - r_i_i_C(2) * t272) * qJD(6)) * t268;
t244 = t270 * t293 - t326 * t272;
t243 = t326 * t270 + t272 * t293;
t242 = t270 * t328 + t288 * t272;
t241 = t288 * t270 - t272 * t328;
t1 = [t244 * r_i_i_C(1) + t243 * r_i_i_C(2) - t330 * t282 + t287 * t284 (-t262 - t311 - t322) * t315 + t286 * t284 + t290, t314 (-r_i_i_C(3) * t316 - t315 * t323) * t267 + (-r_i_i_C(3) * t315 - t284 * t294) * t268 + t290, t268 * t316 - t307, t241 * r_i_i_C(1) + t242 * r_i_i_C(2); -t242 * r_i_i_C(1) + t241 * r_i_i_C(2) + t287 * t282 + t330 * t284, t286 * t282 + (t262 + t291) * t314 + t292, t315, t282 * t289 + t291 * t314 + t292, t268 * t317 + t306, -t243 * r_i_i_C(1) + t244 * r_i_i_C(2); 0, t257 + t285, 0, t285, t277 * t267 (-t272 * t319 + t305) * r_i_i_C(2) - t332 * r_i_i_C(1);];
JaD_transl  = t1;
