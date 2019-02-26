% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR10_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR10_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:43
% EndTime: 2019-02-26 21:33:43
% DurationCPUTime: 0.27s
% Computational Cost: add. (420->60), mult. (519->87), div. (0->0), fcn. (413->10), ass. (0->51)
t242 = sin(qJ(2));
t244 = cos(qJ(2));
t240 = pkin(10) + qJ(5);
t237 = cos(t240);
t282 = pkin(5) * t237;
t255 = qJD(5) * t282 + qJD(3);
t265 = t244 * qJD(4);
t236 = sin(t240);
t275 = qJ(3) + pkin(5) * t236 + sin(pkin(10)) * pkin(4);
t264 = pkin(2) + r_i_i_C(3) + pkin(9) + pkin(8) + qJ(4);
t285 = t264 * t242;
t293 = qJD(2) * (-t275 * t244 + t285) - (pkin(7) + t282 + cos(pkin(10)) * pkin(4) + pkin(3)) * qJD(1) - t255 * t242 - t265;
t238 = qJ(6) + t240;
t234 = sin(t238);
t235 = cos(t238);
t292 = r_i_i_C(1) * t234 + r_i_i_C(2) * t235;
t291 = r_i_i_C(1) * t235 - r_i_i_C(2) * t234;
t253 = t275 + t292;
t289 = -t253 * t244 + t285;
t243 = sin(qJ(1));
t241 = qJD(5) + qJD(6);
t261 = t241 * t242 + qJD(1);
t288 = t243 * t261;
t245 = cos(qJ(1));
t287 = t245 * t261;
t272 = qJD(1) * t242;
t260 = -t241 - t272;
t268 = qJD(2) * t244;
t263 = t243 * t268;
t251 = t260 * t245 - t263;
t225 = t234 * t288 + t251 * t235;
t226 = t251 * t234 - t235 * t288;
t274 = -t225 * r_i_i_C(1) + t226 * r_i_i_C(2);
t267 = qJD(2) * t245;
t262 = t244 * t267;
t250 = t260 * t243 + t262;
t227 = -t234 * t287 + t250 * t235;
t228 = t250 * t234 + t235 * t287;
t273 = t227 * r_i_i_C(1) - t228 * r_i_i_C(2);
t271 = qJD(1) * t243;
t270 = qJD(1) * t245;
t269 = qJD(2) * t242;
t266 = qJD(5) * t236;
t259 = qJD(5) + t272;
t257 = t264 * t244;
t254 = t236 * (-qJD(5) * t242 - qJD(1));
t252 = t241 * t244 * t292 + t291 * t269;
t249 = t241 * t291 + t255;
t247 = -pkin(5) * t266 + (-t275 * t242 - pkin(1) - t257) * qJD(1);
t246 = -qJD(4) * t242 + t249 * t244 + (-t253 * t242 - t257) * qJD(2);
t1 = [t226 * r_i_i_C(1) + t225 * r_i_i_C(2) + t293 * t243 + t247 * t245, t246 * t245 + t271 * t289, -t242 * t271 + t262, -t242 * t267 - t244 * t271 (t245 * t254 + (-t259 * t243 + t262) * t237) * pkin(5) + t273, t273; t228 * r_i_i_C(1) + t227 * r_i_i_C(2) + t247 * t243 - t293 * t245, t246 * t243 - t270 * t289, t242 * t270 + t263, -t243 * t269 + t244 * t270 (t243 * t254 + (t259 * t245 + t263) * t237) * pkin(5) + t274, t274; 0, -qJD(2) * t289 + t249 * t242 + t265, t269, t268 (t237 * t269 + t244 * t266) * pkin(5) + t252, t252;];
JaD_transl  = t1;
