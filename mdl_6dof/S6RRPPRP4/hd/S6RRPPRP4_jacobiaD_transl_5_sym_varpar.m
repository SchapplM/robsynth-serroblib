% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRP4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:49
% EndTime: 2019-02-26 21:26:49
% DurationCPUTime: 0.40s
% Computational Cost: add. (203->66), mult. (659->110), div. (0->0), fcn. (600->8), ass. (0->48)
t241 = sin(qJ(2));
t238 = sin(pkin(9));
t239 = cos(pkin(9));
t240 = sin(qJ(5));
t243 = cos(qJ(5));
t275 = pkin(3) + pkin(4);
t251 = t243 * r_i_i_C(1) - t240 * r_i_i_C(2) + t275;
t252 = t240 * r_i_i_C(1) + t243 * r_i_i_C(2) + qJ(4);
t276 = t252 * t238 + t251 * t239 + pkin(2);
t244 = cos(qJ(2));
t263 = -r_i_i_C(3) - pkin(8) + qJ(3);
t278 = t263 * t244;
t281 = t276 * t241 - t278;
t264 = t241 * qJD(3);
t280 = (-t241 * pkin(2) + t278) * qJD(2) + t264;
t253 = t238 * t240 + t239 * t243;
t254 = t238 * t243 - t239 * t240;
t249 = t254 * r_i_i_C(1) - t253 * r_i_i_C(2);
t277 = t238 * qJD(4) + t249 * qJD(5);
t242 = sin(qJ(1));
t273 = t242 * t244;
t245 = cos(qJ(1));
t272 = t245 * t238;
t271 = t245 * t239;
t270 = qJD(1) * t242;
t269 = qJD(1) * t245;
t268 = qJD(2) * t241;
t267 = qJD(2) * t244;
t266 = qJD(2) * t245;
t262 = t244 * t272;
t261 = t242 * t268;
t260 = t241 * t266;
t259 = t263 * t241;
t230 = t238 * t273 + t271;
t231 = t239 * t273 - t272;
t256 = -t230 * t243 + t231 * t240;
t255 = t230 * t240 + t231 * t243;
t233 = t242 * t238 + t244 * t271;
t250 = -pkin(2) * t244 - pkin(1) - t259;
t246 = t244 * qJD(3) - t277 * t241 + (-t244 * t276 - t259) * qJD(2);
t232 = -t242 * t239 + t262;
t229 = t233 * qJD(1) - t239 * t261;
t228 = qJD(1) * t262 - t238 * t261 - t239 * t270;
t227 = -t231 * qJD(1) - t239 * t260;
t226 = t230 * qJD(1) + t238 * t260;
t225 = -t226 * t240 + t227 * t243 + (t232 * t243 - t233 * t240) * qJD(5);
t224 = -t226 * t243 - t227 * t240 + (-t232 * t240 - t233 * t243) * qJD(5);
t1 = [-t230 * qJD(4) - t252 * t228 - t251 * t229 + (t256 * r_i_i_C(1) + t255 * r_i_i_C(2)) * qJD(5) - t280 * t242 + (-t242 * pkin(7) + t250 * t245) * qJD(1), t246 * t245 + t270 * t281, -t241 * t270 + t244 * t266, -t226, t224 * r_i_i_C(1) - t225 * r_i_i_C(2), 0; t225 * r_i_i_C(1) + t224 * r_i_i_C(2) - t226 * qJ(4) + t232 * qJD(4) + t275 * t227 + t280 * t245 + (pkin(7) * t245 + t250 * t242) * qJD(1), t246 * t242 - t269 * t281, t241 * t269 + t242 * t267, t228 (t228 * t243 - t229 * t240) * r_i_i_C(1) + (-t228 * t240 - t229 * t243) * r_i_i_C(2) + (-t255 * r_i_i_C(1) + t256 * r_i_i_C(2)) * qJD(5), 0; 0, -qJD(2) * t281 + t277 * t244 + t264, t268, t238 * t267 (-t253 * r_i_i_C(1) - t254 * r_i_i_C(2)) * t241 * qJD(5) + t249 * t267, 0;];
JaD_transl  = t1;
