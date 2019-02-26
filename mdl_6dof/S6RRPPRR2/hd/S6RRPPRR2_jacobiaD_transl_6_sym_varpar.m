% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:28:45
% EndTime: 2019-02-26 21:28:46
% DurationCPUTime: 0.32s
% Computational Cost: add. (525->64), mult. (475->91), div. (0->0), fcn. (380->12), ass. (0->56)
t253 = pkin(11) + qJ(5);
t249 = cos(t253);
t242 = pkin(5) * t249 + cos(pkin(11)) * pkin(4) + pkin(3);
t255 = qJ(2) + pkin(10);
t248 = sin(t255);
t250 = cos(t255);
t293 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(4);
t267 = t293 * t250 - sin(qJ(2)) * pkin(2);
t247 = sin(t253);
t292 = pkin(5) * qJD(5);
t279 = t247 * t292;
t280 = t248 * qJD(4);
t298 = pkin(5) * t247;
t307 = (-t242 * t248 + t267) * qJD(2) + (t298 + sin(pkin(11)) * pkin(4) + qJ(3) + pkin(7)) * qJD(1) - t250 * t279 + t280;
t251 = qJ(6) + t253;
t244 = sin(t251);
t296 = r_i_i_C(2) * t244;
t245 = cos(t251);
t297 = r_i_i_C(1) * t245;
t268 = t242 - t296 + t297;
t263 = -t268 * t248 + t267;
t260 = cos(qJ(1));
t254 = qJD(5) + qJD(6);
t274 = t250 * t254 - qJD(1);
t305 = t260 * t274;
t303 = -t293 * t248 - cos(qJ(2)) * pkin(2);
t295 = r_i_i_C(2) * t245;
t271 = r_i_i_C(1) * t244 + t295;
t302 = t271 * t254 + t279;
t286 = qJD(1) * t250;
t273 = -t254 + t286;
t258 = sin(qJ(1));
t282 = qJD(2) * t258;
t277 = t248 * t282;
t301 = t273 * t260 - t277;
t290 = t248 * t254;
t281 = qJD(2) * t260;
t276 = t248 * t281;
t264 = t273 * t258 + t276;
t237 = t264 * t244 - t245 * t305;
t238 = t244 * t305 + t264 * t245;
t289 = t237 * r_i_i_C(1) + t238 * r_i_i_C(2);
t270 = t274 * t258;
t239 = t301 * t244 + t245 * t270;
t240 = t244 * t270 - t301 * t245;
t288 = -t239 * r_i_i_C(1) + t240 * r_i_i_C(2);
t285 = qJD(1) * t258;
t284 = qJD(1) * t260;
t283 = qJD(2) * t250;
t278 = t249 * t292;
t272 = -qJD(5) + t286;
t269 = t249 * (-qJD(5) * t250 + qJD(1));
t262 = t278 + qJD(3) + (-t242 * t250 - pkin(1) + t303) * qJD(1);
t261 = qJD(4) * t250 + t302 * t248 + (-t268 * t250 + t303) * qJD(2);
t241 = t290 * t296;
t1 = [t240 * r_i_i_C(1) + t239 * r_i_i_C(2) - t307 * t258 + t262 * t260, t261 * t260 - t263 * t285, t284, -t248 * t285 + t250 * t281 (t260 * t269 + (t272 * t258 + t276) * t247) * pkin(5) + t289, t289; -t238 * r_i_i_C(1) + t237 * r_i_i_C(2) + t262 * t258 + t307 * t260, t261 * t258 + t263 * t284, t285, t248 * t284 + t250 * t282 (t258 * t269 + (-t272 * t260 + t277) * t247) * pkin(5) + t288, t288; 0, t263 * qJD(2) - t302 * t250 + t280, 0, qJD(2) * t248, t241 + (-t254 * t297 - t278) * t248 + (-t271 - t298) * t283, -t283 * t295 + t241 + (-t244 * t283 - t245 * t290) * r_i_i_C(1);];
JaD_transl  = t1;
