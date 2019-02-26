% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:37
% EndTime: 2019-02-26 20:18:37
% DurationCPUTime: 0.18s
% Computational Cost: add. (411->54), mult. (480->90), div. (0->0), fcn. (433->12), ass. (0->49)
t236 = sin(pkin(12));
t238 = cos(pkin(12));
t241 = sin(qJ(2));
t239 = cos(pkin(6));
t243 = cos(qJ(2));
t259 = t239 * t243;
t271 = -t236 * t241 + t238 * t259;
t235 = qJ(3) + qJ(4);
t230 = sin(t235);
t270 = pkin(4) * t230;
t269 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8);
t268 = pkin(3) * qJD(3);
t232 = qJ(5) + t235;
t227 = sin(t232);
t233 = qJD(3) + qJD(4);
t229 = qJD(5) + t233;
t267 = t227 * t229;
t228 = cos(t232);
t266 = t228 * t229;
t231 = cos(t235);
t265 = t231 * t233;
t237 = sin(pkin(6));
t264 = t236 * t237;
t262 = t237 * t238;
t261 = t237 * t241;
t260 = t239 * t241;
t219 = t236 * t243 + t238 * t260;
t214 = t271 * qJD(2);
t251 = t229 * t262 - t214;
t258 = (-t219 * t266 + t227 * t251) * r_i_i_C(1) + (t219 * t267 + t228 * t251) * r_i_i_C(2);
t247 = t236 * t260 - t238 * t243;
t248 = t236 * t259 + t238 * t241;
t216 = t248 * qJD(2);
t250 = -t229 * t264 + t216;
t257 = (t227 * t250 + t247 * t266) * r_i_i_C(1) + (t228 * t250 - t247 * t267) * r_i_i_C(2);
t255 = qJD(2) * t243;
t252 = t237 * t255;
t246 = -t229 * t239 - t252;
t254 = t229 * t261;
t256 = (t227 * t246 - t228 * t254) * r_i_i_C(1) + (t227 * t254 + t228 * t246) * r_i_i_C(2);
t242 = cos(qJ(3));
t249 = pkin(3) * t242 + pkin(4) * t231 + r_i_i_C(1) * t228 - r_i_i_C(2) * t227 + pkin(2);
t245 = qJD(2) * t249;
t240 = sin(qJ(3));
t222 = -t233 * t270 - t240 * t268;
t244 = t222 + (-r_i_i_C(1) * t227 - r_i_i_C(2) * t228) * t229;
t225 = -pkin(3) * t240 - t270;
t223 = -pkin(4) * t265 - t242 * t268;
t1 = [0, -t216 * t269 - t244 * t248 + t245 * t247, -t216 * t225 + t222 * t264 - t223 * t247 + t257 (t247 * t265 + (-t233 * t264 + t216) * t230) * pkin(4) + t257, t257, 0; 0, t214 * t269 - t219 * t245 + t244 * t271, t214 * t225 + t219 * t223 - t222 * t262 + t258 (-t219 * t265 + (t233 * t262 - t214) * t230) * pkin(4) + t258, t258, 0; 0 (t244 * t243 + (-t241 * t249 + t243 * t269) * qJD(2)) * t237, t239 * t222 + (t223 * t241 + t225 * t255) * t237 + t256 (-t261 * t265 + (-t233 * t239 - t252) * t230) * pkin(4) + t256, t256, 0;];
JaD_transl  = t1;
