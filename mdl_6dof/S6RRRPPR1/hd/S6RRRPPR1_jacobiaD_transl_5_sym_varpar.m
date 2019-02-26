% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:03:29
% EndTime: 2019-02-26 22:03:29
% DurationCPUTime: 0.21s
% Computational Cost: add. (340->50), mult. (317->56), div. (0->0), fcn. (230->10), ass. (0->45)
t236 = qJ(2) + qJ(3);
t231 = pkin(10) + t236;
t230 = cos(t231);
t270 = r_i_i_C(3) + qJ(5);
t257 = t270 * t230;
t229 = sin(t231);
t228 = t229 * qJD(5);
t232 = sin(t236);
t235 = qJD(2) + qJD(3);
t237 = sin(pkin(11));
t238 = cos(pkin(11));
t271 = r_i_i_C(2) * t237;
t280 = r_i_i_C(1) * t238 + pkin(4);
t251 = t280 - t271;
t239 = sin(qJ(2));
t269 = pkin(2) * qJD(2);
t263 = t239 * t269;
t273 = pkin(3) * t235;
t281 = (-t251 * t229 + t257) * t235 + (r_i_i_C(1) * t237 + r_i_i_C(2) * t238 + pkin(7) + pkin(8) + qJ(4)) * qJD(1) - t232 * t273 + t228 - t263;
t268 = t230 * t235;
t279 = qJD(5) * t230 + t268 * t271;
t275 = pkin(3) * t232;
t233 = cos(t236);
t274 = pkin(3) * t233;
t240 = sin(qJ(1));
t267 = qJD(1) * t240;
t242 = cos(qJ(1));
t266 = qJD(1) * t242;
t264 = t229 * t271;
t261 = t229 * t266;
t259 = t229 * t267;
t258 = t270 * t229;
t254 = t279 * t242 + t280 * t259;
t253 = t280 * t229;
t252 = t279 * t240 + t266 * t257 + t261 * t271;
t249 = -t257 - t264;
t248 = -t253 - t275;
t247 = -t230 * t280 - t258;
t241 = cos(qJ(2));
t246 = -t233 * t273 + t247 * t235 - t241 * t269;
t245 = t235 * (t247 - t274);
t244 = t228 + t270 * t268 + (t248 + t264) * t235;
t243 = qJD(4) + (-pkin(2) * t241 - t251 * t230 - pkin(1) - t258 - t274) * qJD(1);
t225 = -pkin(2) * t239 - t275;
t1 = [-t281 * t240 + t243 * t242 (-t225 + t249) * t267 + t246 * t242 + t254 (t249 + t275) * t267 + t242 * t245 + t254, t266, t242 * t268 - t259, 0; t243 * t240 + t281 * t242 (t225 - t253) * t266 + t246 * t240 + t252, t240 * t245 + t248 * t266 + t252, t267, t240 * t268 + t261, 0; 0, t244 - t263, t244, 0, t235 * t229, 0;];
JaD_transl  = t1;
