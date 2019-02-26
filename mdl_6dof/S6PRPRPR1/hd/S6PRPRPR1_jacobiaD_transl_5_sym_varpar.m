% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:46:28
% EndTime: 2019-02-26 19:46:28
% DurationCPUTime: 0.24s
% Computational Cost: add. (201->53), mult. (547->99), div. (0->0), fcn. (552->12), ass. (0->43)
t261 = cos(pkin(6));
t256 = sin(pkin(11));
t259 = cos(pkin(11));
t264 = sin(qJ(2));
t266 = cos(qJ(2));
t270 = t266 * t256 + t264 * t259;
t242 = t270 * t261;
t276 = qJD(2) * t266;
t277 = qJD(2) * t264;
t284 = t256 * t277 - t259 * t276;
t283 = -r_i_i_C(3) - qJ(5) - pkin(8);
t282 = pkin(2) * qJD(2);
t257 = sin(pkin(10));
t258 = sin(pkin(6));
t281 = t257 * t258;
t260 = cos(pkin(10));
t280 = t258 * t260;
t263 = sin(qJ(4));
t279 = t258 * t263;
t278 = t261 * t264;
t237 = t284 * t261;
t244 = -t256 * t276 - t259 * t277;
t273 = t260 * t237 - t257 * t244;
t272 = t257 * t237 + t260 * t244;
t245 = t264 * t256 - t266 * t259;
t232 = t260 * t242 - t257 * t245;
t271 = t257 * t242 + t260 * t245;
t255 = qJ(4) + pkin(12);
t253 = sin(t255);
t254 = cos(t255);
t265 = cos(qJ(4));
t269 = t265 * pkin(4) + t254 * r_i_i_C(1) - t253 * r_i_i_C(2) + pkin(3);
t240 = t270 * t258;
t268 = t263 * pkin(4) + t253 * r_i_i_C(1) + t254 * r_i_i_C(2);
t267 = qJD(4) * t268;
t243 = t245 * qJD(2);
t241 = t245 * t261;
t238 = qJD(2) * t242;
t236 = qJD(2) * t240;
t235 = t284 * t258;
t228 = t257 * t238 + t260 * t243;
t225 = -t260 * t238 + t257 * t243;
t1 = [0, -t271 * qJD(5) - t283 * t272 + (t257 * t278 - t260 * t266) * t282 + t269 * t228 - (t257 * t241 - t260 * t270) * t267, 0, -t268 * t272 + ((-t253 * t281 + t254 * t271) * r_i_i_C(1) + (-t253 * t271 - t254 * t281) * r_i_i_C(2) + (-t257 * t279 + t265 * t271) * pkin(4)) * qJD(4), -t228, 0; 0, t232 * qJD(5) + t283 * t273 + (-t257 * t266 - t260 * t278) * t282 + t269 * t225 - (-t260 * t241 - t257 * t270) * t267, 0, t268 * t273 + ((-t232 * t254 + t253 * t280) * r_i_i_C(1) + (t232 * t253 + t254 * t280) * r_i_i_C(2) + (-t232 * t265 + t260 * t279) * pkin(4)) * qJD(4), -t225, 0; 0, t240 * qJD(5) + t283 * t235 - t269 * t236 + (-pkin(2) * t277 + t245 * t267) * t258, 0, t268 * t235 + ((-t240 * t254 - t253 * t261) * r_i_i_C(1) + (t240 * t253 - t254 * t261) * r_i_i_C(2) + (-t240 * t265 - t261 * t263) * pkin(4)) * qJD(4), t236, 0;];
JaD_transl  = t1;
