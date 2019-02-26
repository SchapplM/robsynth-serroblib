% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:35:25
% EndTime: 2019-02-26 20:35:26
% DurationCPUTime: 0.28s
% Computational Cost: add. (491->62), mult. (430->96), div. (0->0), fcn. (340->11), ass. (0->51)
t263 = cos(qJ(5));
t250 = pkin(5) * t263 + pkin(4);
t257 = pkin(11) + qJ(4);
t251 = sin(t257);
t253 = cos(t257);
t293 = r_i_i_C(3) + pkin(9) + pkin(8);
t279 = t293 * t253;
t262 = sin(qJ(5));
t292 = pkin(5) * qJD(5);
t282 = t262 * t292;
t301 = (-t250 * t251 + t279) * qJD(4) - t253 * t282;
t260 = qJ(5) + qJ(6);
t256 = cos(t260);
t258 = qJD(5) + qJD(6);
t286 = qJD(1) * t253;
t276 = -t258 + t286;
t299 = t256 * t276;
t298 = t262 * (-qJD(5) + t286);
t277 = t253 * t258 - qJD(1);
t255 = sin(t260);
t284 = qJD(4) * t255;
t297 = -t251 * t284 + t277 * t256;
t272 = r_i_i_C(1) * t255 + r_i_i_C(2) * t256;
t296 = -t272 * t258 - t282;
t295 = pkin(5) * t262;
t294 = r_i_i_C(2) * t255;
t290 = t256 * t258;
t259 = qJ(1) + pkin(10);
t252 = sin(t259);
t254 = cos(t259);
t271 = t276 * t255;
t244 = t252 * t271 - t297 * t254;
t283 = qJD(4) * t256;
t267 = t251 * t283 + t277 * t255;
t245 = t252 * t299 + t267 * t254;
t289 = t244 * r_i_i_C(1) + t245 * r_i_i_C(2);
t246 = t297 * t252 + t254 * t271;
t247 = t267 * t252 - t254 * t299;
t288 = -t246 * r_i_i_C(1) + t247 * r_i_i_C(2);
t287 = qJD(1) * t252;
t285 = qJD(1) * t254;
t281 = t263 * t292;
t278 = pkin(7) + qJ(3) + t295;
t273 = qJD(3) + t281;
t270 = r_i_i_C(1) * t256 + t250 - t294;
t269 = qJD(4) * t270;
t268 = -t250 * t253 - t293 * t251 - cos(pkin(11)) * pkin(3) - pkin(2);
t266 = qJD(4) * t251 * t262 + (-qJD(5) * t253 + qJD(1)) * t263;
t265 = -t293 * qJD(4) - t296;
t248 = t251 * t258 * t294;
t1 = [t247 * r_i_i_C(1) + t246 * r_i_i_C(2) + t273 * t254 - t301 * t252 + (-cos(qJ(1)) * pkin(1) - t278 * t252 + t268 * t254) * qJD(1), 0, t285 (-t254 * t269 - t293 * t287) * t253 + (t265 * t254 + t270 * t287) * t251 (t252 * t298 + t266 * t254) * pkin(5) + t289, t289; -t245 * r_i_i_C(1) + t244 * r_i_i_C(2) + t273 * t252 + t301 * t254 + (-sin(qJ(1)) * pkin(1) + t278 * t254 + t268 * t252) * qJD(1), 0, t287 (-t252 * t269 + t293 * t285) * t253 + (t265 * t252 - t270 * t285) * t251 (t266 * t252 - t254 * t298) * pkin(5) + t288, t288; 0, 0, 0, t296 * t253 + (-t270 * t251 + t279) * qJD(4), t248 + (-r_i_i_C(1) * t290 - t281) * t251 + (-t272 - t295) * t253 * qJD(4), -r_i_i_C(2) * t253 * t283 + t248 + (-t251 * t290 - t253 * t284) * r_i_i_C(1);];
JaD_transl  = t1;
