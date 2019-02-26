% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPPR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:23:46
% EndTime: 2019-02-26 21:23:46
% DurationCPUTime: 0.35s
% Computational Cost: add. (222->69), mult. (707->109), div. (0->0), fcn. (637->8), ass. (0->49)
t247 = cos(qJ(2));
t241 = sin(pkin(9));
t242 = cos(pkin(9));
t243 = sin(qJ(6));
t246 = cos(qJ(6));
t280 = pkin(4) + pkin(5);
t255 = t246 * r_i_i_C(1) - t243 * r_i_i_C(2) + t280;
t256 = t243 * r_i_i_C(1) + t246 * r_i_i_C(2) + qJ(5);
t282 = t255 * t241 - t256 * t242 + qJ(3);
t244 = sin(qJ(2));
t265 = r_i_i_C(3) + pkin(8) - qJ(4) - pkin(2);
t285 = t265 * t244;
t250 = t282 * t247 + t285;
t269 = t247 * qJD(4);
t287 = (t247 * qJ(3) + t285) * qJD(2) + (pkin(3) + pkin(7)) * qJD(1) + t244 * qJD(3) + t269;
t245 = sin(qJ(1));
t271 = qJD(2) * t247;
t267 = t245 * t271;
t248 = cos(qJ(1));
t273 = qJD(1) * t248;
t284 = t244 * t273 + t267;
t279 = t245 * t241;
t278 = t245 * t242;
t276 = t248 * t241;
t275 = t248 * t242;
t274 = qJD(1) * t245;
t272 = qJD(2) * t244;
t270 = qJD(2) * t248;
t266 = t247 * t270;
t262 = t265 * t247;
t236 = t244 * t278 + t276;
t237 = -t244 * t279 + t275;
t261 = t236 * t246 - t237 * t243;
t260 = t236 * t243 + t237 * t246;
t259 = t241 * t246 - t242 * t243;
t258 = -t241 * t243 - t242 * t246;
t235 = t244 * t276 + t278;
t254 = qJD(1) * (-qJ(3) * t244 - pkin(1) + t262);
t253 = t258 * r_i_i_C(1) - t259 * r_i_i_C(2);
t251 = -t242 * qJD(5) + t253 * qJD(6) + qJD(3);
t249 = -t244 * qJD(4) + t251 * t247 + (-t244 * t282 + t262) * qJD(2);
t234 = -t244 * t275 + t279;
t233 = t237 * qJD(1) + t241 * t266;
t232 = t236 * qJD(1) - t242 * t266;
t231 = t235 * qJD(1) + t241 * t267;
t230 = t241 * t274 - t284 * t242;
t229 = t232 * t243 + t233 * t246 + (t234 * t246 - t235 * t243) * qJD(6);
t228 = t232 * t246 - t233 * t243 + (-t234 * t243 - t235 * t246) * qJD(6);
t1 = [t236 * qJD(5) - t256 * t230 - t255 * t231 + (t261 * r_i_i_C(1) - t260 * r_i_i_C(2)) * qJD(6) + t248 * t254 - t287 * t245, t249 * t248 - t250 * t274, -t244 * t274 + t266, -t244 * t270 - t247 * t274, t232, t228 * r_i_i_C(1) - t229 * r_i_i_C(2); t229 * r_i_i_C(1) + t228 * r_i_i_C(2) + t232 * qJ(5) + t234 * qJD(5) + t280 * t233 + t245 * t254 + t287 * t248, t249 * t245 + t250 * t273, t284, -t245 * t272 + t247 * t273, t230 (t230 * t246 - t231 * t243) * r_i_i_C(1) + (-t230 * t243 - t231 * t246) * r_i_i_C(2) + (t260 * r_i_i_C(1) + t261 * r_i_i_C(2)) * qJD(6); 0, t250 * qJD(2) + t251 * t244 + t269, t272, t271, -t242 * t272 (t259 * r_i_i_C(1) + t258 * r_i_i_C(2)) * t247 * qJD(6) + t253 * t272;];
JaD_transl  = t1;
