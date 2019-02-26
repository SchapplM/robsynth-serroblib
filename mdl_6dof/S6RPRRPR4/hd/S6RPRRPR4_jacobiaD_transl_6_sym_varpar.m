% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR4
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
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:38
% EndTime: 2019-02-26 21:02:38
% DurationCPUTime: 0.31s
% Computational Cost: add. (532->72), mult. (443->102), div. (0->0), fcn. (344->11), ass. (0->57)
t271 = pkin(11) + qJ(6);
t265 = sin(t271);
t267 = cos(t271);
t272 = pkin(10) + qJ(3);
t269 = qJ(4) + t272;
t262 = sin(t269);
t305 = qJD(6) * t262;
t263 = cos(t269);
t273 = qJD(3) + qJD(4);
t312 = t263 * t273;
t327 = t265 * t312 + t267 * t305;
t264 = cos(pkin(11)) * pkin(5) + pkin(4);
t326 = r_i_i_C(1) * t267 + t264;
t261 = t262 * qJD(5);
t266 = sin(t272);
t314 = pkin(3) * qJD(3);
t303 = t266 * t314;
t275 = -pkin(9) - qJ(5);
t315 = r_i_i_C(3) - t275;
t325 = (-t262 * t264 + t315 * t263) * t273 + (pkin(5) * sin(pkin(11)) + pkin(7) + pkin(8) + qJ(2)) * qJD(1) + t261 - t303;
t297 = t265 * t305;
t324 = r_i_i_C(1) * t297 + t327 * r_i_i_C(2) + qJD(5) * t263;
t277 = cos(qJ(1));
t290 = qJD(6) * t263 - qJD(1);
t323 = t277 * t290;
t289 = qJD(1) * t263 - qJD(6);
t276 = sin(qJ(1));
t310 = t273 * t276;
t302 = t262 * t310;
t321 = t289 * t277 - t302;
t319 = pkin(3) * t266;
t317 = r_i_i_C(2) * t265;
t316 = r_i_i_C(3) * t263;
t311 = t263 * t275;
t309 = t273 * t277;
t308 = qJD(1) * t276;
t307 = qJD(1) * t277;
t304 = t262 * t317;
t301 = t262 * t309;
t299 = t262 * t308;
t298 = t262 * t307;
t288 = t326 * t273;
t287 = t290 * t276;
t285 = t275 * t302 + t324 * t276 + t298 * t317 + t307 * t316;
t284 = -t262 * t326 - t311;
t283 = t275 * t301 + t324 * t277 + t326 * t299 + t308 * t311;
t282 = (-r_i_i_C(3) * t262 - t263 * t326) * t273;
t281 = t289 * t276 + t301;
t268 = cos(t272);
t280 = qJD(2) + (-t315 * t262 - t263 * t264 - pkin(3) * t268 - cos(pkin(10)) * pkin(2) - pkin(1)) * qJD(1);
t279 = -t268 * t314 + t282;
t278 = t273 * t304 + r_i_i_C(3) * t312 + t261 - t262 * t288 + (-t273 * t275 + (-r_i_i_C(1) * t265 - r_i_i_C(2) * t267) * qJD(6)) * t263;
t242 = t265 * t287 - t321 * t267;
t241 = t321 * t265 + t267 * t287;
t240 = t265 * t323 + t281 * t267;
t239 = t281 * t265 - t267 * t323;
t1 = [t242 * r_i_i_C(1) + t241 * r_i_i_C(2) - t325 * t276 + t280 * t277, t307 (-t304 - t316 + t319) * t308 + t279 * t277 + t283 (-r_i_i_C(3) * t309 - t308 * t317) * t262 + (-r_i_i_C(3) * t308 - t277 * t288) * t263 + t283, t263 * t309 - t299, t239 * r_i_i_C(1) + t240 * r_i_i_C(2); -t240 * r_i_i_C(1) + t239 * r_i_i_C(2) + t280 * t276 + t325 * t277, t308, t279 * t276 + (t284 - t319) * t307 + t285, t276 * t282 + t284 * t307 + t285, t263 * t310 + t298, -t241 * r_i_i_C(1) + t242 * r_i_i_C(2); 0, 0, t278 - t303, t278, t273 * t262 (-t267 * t312 + t297) * r_i_i_C(2) - t327 * r_i_i_C(1);];
JaD_transl  = t1;
