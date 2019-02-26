% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR11_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR11_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:20
% EndTime: 2019-02-26 21:34:20
% DurationCPUTime: 0.19s
% Computational Cost: add. (268->63), mult. (628->102), div. (0->0), fcn. (595->10), ass. (0->44)
t251 = pkin(11) + qJ(5);
t249 = sin(t251);
t250 = cos(t251);
t263 = r_i_i_C(1) * t250 - r_i_i_C(2) * t249;
t260 = t263 * qJD(5) + qJD(3);
t283 = pkin(8) + cos(pkin(11)) * pkin(4) + pkin(3);
t253 = sin(pkin(6));
t257 = sin(qJ(1));
t282 = t253 * t257;
t258 = cos(qJ(2));
t281 = t253 * t258;
t259 = cos(qJ(1));
t280 = t253 * t259;
t256 = sin(qJ(2));
t279 = t257 * t256;
t278 = t257 * t258;
t277 = t259 * t256;
t276 = t259 * t258;
t275 = qJD(1) * t257;
t274 = qJD(1) * t259;
t273 = qJD(2) * t256;
t272 = qJD(2) * t258;
t271 = -r_i_i_C(3) - pkin(9) - qJ(4) - pkin(2);
t254 = cos(pkin(6));
t270 = t254 * t279;
t269 = t254 * t276;
t268 = t253 * t275;
t267 = t253 * t274;
t266 = t253 * t273;
t265 = -sin(pkin(11)) * pkin(4) - qJ(3);
t264 = qJD(2) * t254 + qJD(1);
t262 = -t249 * r_i_i_C(1) - t250 * r_i_i_C(2);
t240 = t254 * t278 + t277;
t239 = t254 * t277 + t278;
t261 = -t262 - t265;
t241 = -t270 + t276;
t238 = -t269 + t279;
t237 = -qJD(1) * t270 - t257 * t273 + t264 * t276;
t236 = t240 * qJD(1) + t239 * qJD(2);
t235 = t239 * qJD(1) + t240 * qJD(2);
t234 = -qJD(1) * t269 - t259 * t272 + t264 * t279;
t233 = t250 * t267 - t234 * t249 + (t240 * t250 - t249 * t282) * qJD(5);
t232 = -t249 * t267 - t234 * t250 + (-t240 * t249 - t250 * t282) * qJD(5);
t1 = [-pkin(1) * t274 - t239 * qJD(4) - t260 * t238 - t261 * t236 + t271 * t237 + (t262 * t259 * qJD(5) + (-t263 - t283) * t275) * t253, -t240 * qJD(4) - t271 * t234 - t261 * t235 + t260 * t241, -t234, -t235, t232 * r_i_i_C(1) - t233 * r_i_i_C(2), 0; t233 * r_i_i_C(1) + t232 * r_i_i_C(2) + t240 * qJD(3) + t241 * qJD(4) + t265 * t234 + t271 * t235 + (-pkin(1) * t257 + t283 * t280) * qJD(1), -t238 * qJD(4) + t271 * t236 + t261 * t237 + t260 * t239, t236, t237 (t236 * t250 - t249 * t268) * r_i_i_C(1) + (-t236 * t249 - t250 * t268) * r_i_i_C(2) + ((-t238 * t249 + t250 * t280) * r_i_i_C(1) + (-t238 * t250 - t249 * t280) * r_i_i_C(2)) * qJD(5), 0; 0 (qJD(4) * t258 + t260 * t256 + (t271 * t256 + t261 * t258) * qJD(2)) * t253, t266, t253 * t272, t263 * t266 + ((t249 * t281 - t250 * t254) * r_i_i_C(1) + (t249 * t254 + t250 * t281) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
