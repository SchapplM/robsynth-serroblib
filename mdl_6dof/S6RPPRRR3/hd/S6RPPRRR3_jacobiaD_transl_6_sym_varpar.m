% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:35:58
% EndTime: 2019-02-26 20:35:58
% DurationCPUTime: 0.26s
% Computational Cost: add. (391->53), mult. (432->83), div. (0->0), fcn. (340->10), ass. (0->48)
t238 = sin(qJ(4));
t240 = cos(qJ(4));
t234 = qJD(5) + qJD(6);
t237 = sin(qJ(5));
t274 = pkin(5) * t237;
t261 = qJD(5) * t274;
t236 = qJ(5) + qJ(6);
t232 = sin(t236);
t233 = cos(t236);
t285 = r_i_i_C(1) * t232 + r_i_i_C(2) * t233;
t244 = t285 * t234 + t261;
t239 = cos(qJ(5));
t229 = t239 * pkin(5) + pkin(4);
t284 = -r_i_i_C(1) * t233 + r_i_i_C(2) * t232;
t248 = t229 - t284;
t269 = r_i_i_C(3) + pkin(9) + pkin(8);
t276 = t269 * t240;
t287 = (-t248 * t238 + t276) * qJD(4) - t244 * t240;
t257 = t269 * t238;
t286 = t248 * t240 + t257;
t283 = -t229 * t238 - qJ(3) + t276;
t265 = qJD(1) * t238;
t253 = t234 + t265;
t281 = t232 * t253;
t280 = t233 * t253;
t263 = qJD(4) * t240;
t275 = (qJD(5) * t238 + qJD(1)) * t239 + t237 * t263;
t235 = qJ(1) + pkin(10);
t230 = sin(t235);
t231 = cos(t235);
t254 = -t234 * t238 - qJD(1);
t246 = -t232 * t263 + t254 * t233;
t222 = t230 * t281 + t246 * t231;
t245 = t254 * t232 + t233 * t263;
t223 = -t230 * t280 + t245 * t231;
t267 = -t222 * r_i_i_C(1) + t223 * r_i_i_C(2);
t224 = t246 * t230 - t231 * t281;
t225 = t245 * t230 + t231 * t280;
t266 = t224 * r_i_i_C(1) - t225 * r_i_i_C(2);
t264 = qJD(4) * t238;
t262 = qJD(5) * t239;
t260 = pkin(5) * t262;
t252 = -pkin(2) - pkin(7) - t274;
t249 = (-qJD(5) - t265) * t237;
t247 = t284 * t234 * t240 + t285 * t264;
t243 = qJD(1) * t286;
t242 = -t238 * t261 + qJD(3) + (t229 * t240 + t257) * qJD(4);
t1 = [-t230 * t260 + t223 * r_i_i_C(1) + t222 * r_i_i_C(2) + t242 * t231 + (-cos(qJ(1)) * pkin(1) + t252 * t231 + t283 * t230) * qJD(1), 0, qJD(1) * t231, t287 * t230 + t231 * t243 (-t275 * t230 + t231 * t249) * pkin(5) + t266, t266; t231 * t260 + t225 * r_i_i_C(1) + t224 * r_i_i_C(2) + t242 * t230 + (-sin(qJ(1)) * pkin(1) + t252 * t230 - t283 * t231) * qJD(1), 0, qJD(1) * t230, t230 * t243 - t287 * t231 (t230 * t249 + t275 * t231) * pkin(5) + t267, t267; 0, 0, 0, -t286 * qJD(4) + t244 * t238 (t237 * t264 - t240 * t262) * pkin(5) + t247, t247;];
JaD_transl  = t1;
