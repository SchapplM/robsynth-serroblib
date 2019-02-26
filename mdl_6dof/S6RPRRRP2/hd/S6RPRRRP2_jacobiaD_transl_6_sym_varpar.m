% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:08:38
% EndTime: 2019-02-26 21:08:38
% DurationCPUTime: 0.28s
% Computational Cost: add. (516->69), mult. (526->98), div. (0->0), fcn. (412->10), ass. (0->57)
t246 = qJ(1) + pkin(10);
t239 = sin(t246);
t247 = qJ(4) + qJ(5);
t241 = sin(t247);
t242 = cos(t247);
t245 = qJD(4) + qJD(5);
t251 = cos(qJ(3));
t265 = t245 * t251 - qJD(1);
t249 = sin(qJ(3));
t272 = qJD(3) * t249;
t287 = -t241 * t272 + t242 * t265;
t294 = t287 * t239;
t250 = cos(qJ(4));
t237 = t250 * pkin(4) + pkin(5) * t242;
t235 = pkin(3) + t237;
t284 = r_i_i_C(2) * t241;
t258 = r_i_i_C(1) * t242 + t235 - t284;
t281 = r_i_i_C(3) + qJ(6) + pkin(9) + pkin(8);
t289 = t281 * t251;
t253 = -t249 * t258 + t289;
t293 = qJD(1) * t253;
t248 = sin(qJ(4));
t280 = pkin(4) * qJD(4);
t285 = pkin(5) * t241;
t232 = -t245 * t285 - t248 * t280;
t270 = t249 * qJD(6);
t292 = (-t235 * t249 + t289) * qJD(3) + t251 * t232 + t270;
t283 = r_i_i_C(2) * t242;
t263 = r_i_i_C(1) * t241 + t283;
t288 = t245 * t263 - t232;
t286 = -pkin(5) - r_i_i_C(1);
t236 = t248 * pkin(4) + t285;
t282 = pkin(7) + t236;
t278 = t242 * t245;
t277 = t245 * t249;
t240 = cos(t246);
t273 = qJD(1) * t251;
t264 = -t245 + t273;
t260 = t264 * t241;
t228 = t239 * t260 - t240 * t287;
t254 = t241 * t265 + t242 * t272;
t229 = t239 * t242 * t264 + t240 * t254;
t276 = t228 * r_i_i_C(1) + t229 * r_i_i_C(2);
t230 = t240 * t260 + t294;
t261 = t240 * t264;
t231 = t239 * t254 - t242 * t261;
t275 = -t230 * r_i_i_C(1) + t231 * r_i_i_C(2);
t274 = qJD(1) * t249;
t271 = qJD(3) * t251;
t268 = t281 * t249;
t262 = t236 * t273 + t232;
t257 = -t235 * t251 - pkin(2) - t268;
t233 = pkin(5) * t278 + t250 * t280;
t256 = qJD(1) * t237 - t233 * t251 + t236 * t272;
t252 = qJD(6) * t251 + t288 * t249 + (-t251 * t258 - t268) * qJD(3);
t234 = t277 * t284;
t1 = [t231 * r_i_i_C(1) + t230 * r_i_i_C(2) + t240 * t233 - t292 * t239 + (-cos(qJ(1)) * pkin(1) - t282 * t239 + t257 * t240) * qJD(1), 0, -t239 * t293 + t252 * t240, t239 * t262 + t240 * t256 + t276, pkin(5) * t228 + t276, -t239 * t274 + t240 * t271; -t229 * r_i_i_C(1) + t228 * r_i_i_C(2) + t239 * t233 + t292 * t240 + (-sin(qJ(1)) * pkin(1) + t282 * t240 + t257 * t239) * qJD(1), 0, t252 * t239 + t240 * t293, t239 * t256 - t240 * t262 + t275 (-t241 * t261 - t294) * pkin(5) + t275, t239 * t271 + t240 * t274; 0, 0, t253 * qJD(3) - t288 * t251 + t270, t234 + (-r_i_i_C(1) * t278 - t233) * t249 + (-t236 - t263) * t271, t234 + t286 * t242 * t277 + (t241 * t286 - t283) * t271, t272;];
JaD_transl  = t1;
