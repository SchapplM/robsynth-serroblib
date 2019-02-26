% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR11_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR11_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:37
% EndTime: 2019-02-26 21:43:37
% DurationCPUTime: 0.27s
% Computational Cost: add. (439->61), mult. (550->84), div. (0->0), fcn. (432->10), ass. (0->52)
t246 = qJ(4) + pkin(10);
t236 = pkin(5) * cos(t246) + cos(qJ(4)) * pkin(4);
t248 = sin(qJ(2));
t251 = cos(qJ(2));
t230 = t236 * qJD(4);
t271 = qJD(3) + t230;
t272 = t251 * qJD(5);
t235 = sin(qJ(4)) * pkin(4) + pkin(5) * sin(t246);
t281 = qJ(3) + t235;
t270 = pkin(2) + r_i_i_C(3) + pkin(9) + qJ(5) + pkin(8);
t289 = t270 * t248;
t297 = (-t251 * t281 + t289) * qJD(2) - (pkin(7) + pkin(3) + t236) * qJD(1) - t271 * t248 - t272;
t242 = qJ(6) + t246;
t238 = sin(t242);
t239 = cos(t242);
t296 = r_i_i_C(1) * t238 + r_i_i_C(2) * t239;
t260 = t281 + t296;
t295 = -t260 * t251 + t289;
t249 = sin(qJ(1));
t245 = qJD(4) + qJD(6);
t266 = t245 * t248 + qJD(1);
t293 = t249 * t266;
t252 = cos(qJ(1));
t291 = t252 * t266;
t286 = r_i_i_C(1) * t239;
t285 = r_i_i_C(2) * t238;
t278 = qJD(1) * t248;
t265 = -t245 - t278;
t274 = qJD(2) * t251;
t268 = t249 * t274;
t257 = t252 * t265 - t268;
t225 = t238 * t293 + t257 * t239;
t226 = t238 * t257 - t239 * t293;
t280 = -t225 * r_i_i_C(1) + t226 * r_i_i_C(2);
t273 = qJD(2) * t252;
t267 = t251 * t273;
t256 = t249 * t265 + t267;
t227 = -t238 * t291 + t239 * t256;
t228 = t238 * t256 + t239 * t291;
t279 = t227 * r_i_i_C(1) - t228 * r_i_i_C(2);
t277 = qJD(1) * t249;
t276 = qJD(1) * t252;
t275 = qJD(2) * t248;
t269 = t296 * t245 * t251 + t275 * t286;
t263 = t270 * t251;
t261 = t236 * t278 + t230;
t259 = (-t285 + t286) * t245 + t271;
t229 = t235 * qJD(4);
t258 = -qJD(1) * t235 - t229 * t248 + t236 * t274;
t255 = -t229 + (-t248 * t281 - pkin(1) - t263) * qJD(1);
t253 = -qJD(5) * t248 + t259 * t251 + (-t260 * t248 - t263) * qJD(2);
t1 = [t226 * r_i_i_C(1) + t225 * r_i_i_C(2) + t297 * t249 + t255 * t252, t253 * t252 + t295 * t277, -t248 * t277 + t267, -t249 * t261 + t252 * t258 + t279, -t248 * t273 - t251 * t277, t279; t228 * r_i_i_C(1) + t227 * r_i_i_C(2) + t255 * t249 - t297 * t252, t249 * t253 - t276 * t295, t248 * t276 + t268, t249 * t258 + t252 * t261 + t280, -t249 * t275 + t251 * t276, t280; 0, -qJD(2) * t295 + t248 * t259 + t272, t275, t251 * t229 + (t236 - t285) * t275 + t269, t274, -t275 * t285 + t269;];
JaD_transl  = t1;
