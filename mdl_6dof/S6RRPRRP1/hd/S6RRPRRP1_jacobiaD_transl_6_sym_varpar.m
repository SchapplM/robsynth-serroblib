% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP1
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

function JaD_transl = S6RRPRRP1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:46:01
% EndTime: 2019-02-26 21:46:02
% DurationCPUTime: 0.35s
% Computational Cost: add. (520->75), mult. (533->105), div. (0->0), fcn. (407->10), ass. (0->55)
t329 = pkin(5) + r_i_i_C(1);
t272 = qJD(2) + qJD(4);
t275 = sin(qJ(5));
t323 = r_i_i_C(2) * t275;
t336 = t272 * t323 + qJD(6);
t273 = qJ(2) + pkin(10);
t270 = qJ(4) + t273;
t265 = sin(t270);
t266 = cos(t270);
t278 = cos(qJ(5));
t312 = qJD(5) * t278;
t313 = qJD(5) * t275;
t335 = t266 * t336 + (r_i_i_C(2) * t312 + t313 * t329) * t265;
t277 = sin(qJ(1));
t280 = cos(qJ(1));
t297 = qJD(1) * t266 - qJD(5);
t290 = t297 * t280;
t298 = qJD(5) * t266 - qJD(1);
t292 = t298 * t278;
t317 = t272 * t277;
t308 = t265 * t317;
t239 = t277 * t292 + (t290 - t308) * t275;
t267 = t278 * pkin(5) + pkin(4);
t334 = r_i_i_C(1) * t278 + t267;
t260 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t273);
t254 = t260 * qJD(2);
t274 = -qJ(6) - pkin(9);
t321 = r_i_i_C(3) - t274;
t330 = (-pkin(5) * t313 + t321 * t272) * t266 + (pkin(5) * t275 + pkin(7) + pkin(8) + qJ(3)) * qJD(1) - (t267 * t272 - qJD(6)) * t265 + t254;
t322 = r_i_i_C(3) * t266;
t320 = t266 * t272;
t319 = t266 * t274;
t318 = t272 * t265;
t316 = t272 * t280;
t315 = qJD(1) * t277;
t314 = qJD(1) * t280;
t307 = t265 * t316;
t306 = t265 * t315;
t305 = t265 * t314;
t294 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t273);
t293 = t334 * t272;
t291 = t298 * t275;
t289 = -r_i_i_C(2) * t278 - t275 * t329;
t288 = -t265 * t334 - t319;
t287 = t274 * t308 + t335 * t277 + t305 * t323 + t314 * t322;
t286 = t274 * t307 + t335 * t280 + t334 * t306 + t315 * t319;
t285 = (-r_i_i_C(3) * t265 - t266 * t334) * t272;
t284 = t297 * t277 + t307;
t283 = t294 * qJD(2) + t285;
t282 = pkin(5) * t312 + qJD(3) + (-t321 * t265 - t266 * t267 - pkin(1) + t294) * qJD(1);
t237 = t284 * t275 - t280 * t292;
t281 = r_i_i_C(3) * t320 + (t289 * qJD(5) - t272 * t274) * t266 + (-t293 + t336) * t265;
t240 = -t278 * t290 + (t278 * t318 + t291) * t277;
t238 = t284 * t278 + t280 * t291;
t1 = [t240 * r_i_i_C(1) + t239 * r_i_i_C(2) - t330 * t277 + t282 * t280 (-t265 * t323 - t260 - t322) * t315 + t283 * t280 + t286, t314 (-r_i_i_C(3) * t316 - t315 * t323) * t265 + (-r_i_i_C(3) * t315 - t280 * t293) * t266 + t286, t238 * r_i_i_C(2) + t329 * t237, t266 * t316 - t306; -t238 * r_i_i_C(1) + t237 * r_i_i_C(2) + t282 * t277 + t330 * t280, t283 * t277 + (t260 + t288) * t314 + t287, t315, t277 * t285 + t288 * t314 + t287, t240 * r_i_i_C(2) - t239 * t329, t266 * t317 + t305; 0, t254 + t281, 0, t281, t289 * t320 + (-t278 * t329 + t323) * t265 * qJD(5), t318;];
JaD_transl  = t1;
