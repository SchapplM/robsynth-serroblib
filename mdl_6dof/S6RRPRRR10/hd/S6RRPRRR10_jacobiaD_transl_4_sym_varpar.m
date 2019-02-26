% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR10
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR10_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR10_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_jacobiaD_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:18
% EndTime: 2019-02-26 21:59:18
% DurationCPUTime: 0.22s
% Computational Cost: add. (222->61), mult. (494->100), div. (0->0), fcn. (467->10), ass. (0->45)
t258 = sin(qJ(1));
t255 = cos(pkin(6));
t266 = qJD(2) * t255 + qJD(1);
t257 = sin(qJ(2));
t279 = t258 * t257;
t271 = t255 * t279;
t274 = qJD(2) * t257;
t259 = cos(qJ(2));
t260 = cos(qJ(1));
t276 = t260 * t259;
t238 = -qJD(1) * t271 - t258 * t274 + t266 * t276;
t254 = sin(pkin(6));
t280 = t254 * t260;
t267 = qJD(4) * t280;
t284 = t238 - t267;
t252 = pkin(12) + qJ(4);
t250 = sin(t252);
t251 = cos(t252);
t265 = r_i_i_C(1) * t250 + r_i_i_C(2) * t251;
t262 = qJD(4) * t265;
t283 = r_i_i_C(3) + pkin(9) + qJ(3);
t282 = t254 * t257;
t281 = t254 * t258;
t278 = t258 * t259;
t277 = t260 * t257;
t275 = qJD(1) * t254;
t273 = qJD(2) * t259;
t240 = t255 * t277 + t278;
t272 = qJD(4) * t240;
t270 = t255 * t276;
t269 = pkin(3) * sin(pkin(12)) + pkin(8);
t268 = t260 * t275;
t249 = cos(pkin(12)) * pkin(3) + pkin(2);
t264 = t251 * r_i_i_C(1) - t250 * r_i_i_C(2) + t249;
t263 = t255 * t278 + t277;
t261 = t258 * t275 - t272;
t243 = t251 * t267;
t242 = -t271 + t276;
t239 = -t270 + t279;
t237 = t263 * qJD(1) + t240 * qJD(2);
t236 = t240 * qJD(1) + t263 * qJD(2);
t235 = -qJD(1) * t270 - t260 * t273 + t266 * t279;
t234 = t250 * t268 - t236 * t251 + (-t242 * t250 + t251 * t281) * qJD(4);
t233 = t251 * t268 + t236 * t250 + (-t242 * t251 - t250 * t281) * qJD(4);
t1 = [(-t238 * t251 + t250 * t272 + t243) * r_i_i_C(1) + (t284 * t250 + t251 * t272) * r_i_i_C(2) - t238 * t249 - t239 * qJD(3) - t283 * t237 + (-t260 * pkin(1) + (-t265 - t269) * t281) * qJD(1), t242 * qJD(3) + t264 * t235 - t283 * t236 + t262 * t263, -t235, t233 * r_i_i_C(1) - t234 * r_i_i_C(2), 0, 0; t234 * r_i_i_C(1) + t233 * r_i_i_C(2) + t263 * qJD(3) - t236 * t249 - t283 * t235 + (-pkin(1) * t258 + t269 * t280) * qJD(1), t240 * qJD(3) - t264 * t237 + t283 * t238 + t239 * t262, t237, t243 * r_i_i_C(2) + (t261 * r_i_i_C(1) - t238 * r_i_i_C(2)) * t251 + (-t284 * r_i_i_C(1) - t261 * r_i_i_C(2)) * t250, 0, 0; 0 (qJD(3) * t257 - t259 * t262 + (-t264 * t257 + t283 * t259) * qJD(2)) * t254, t254 * t274, -t265 * t254 * t273 + ((-t250 * t255 - t251 * t282) * r_i_i_C(1) + (t250 * t282 - t251 * t255) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
