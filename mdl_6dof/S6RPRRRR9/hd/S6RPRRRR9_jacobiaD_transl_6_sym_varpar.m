% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR9_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR9_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:15
% EndTime: 2019-02-26 21:19:15
% DurationCPUTime: 0.30s
% Computational Cost: add. (582->69), mult. (580->99), div. (0->0), fcn. (458->10), ass. (0->61)
t253 = cos(qJ(4));
t249 = qJ(4) + qJ(5);
t244 = cos(t249);
t247 = qJD(4) + qJD(5);
t284 = t244 * t247;
t285 = pkin(4) * qJD(4);
t232 = pkin(5) * t284 + t253 * t285;
t238 = t253 * pkin(4) + pkin(5) * t244;
t234 = pkin(3) + t238;
t251 = sin(qJ(3));
t254 = cos(qJ(3));
t286 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8);
t270 = t286 * t254;
t299 = (-t234 * t251 - qJ(2) + t270) * qJD(1) - t232;
t245 = qJ(6) + t249;
t240 = sin(t245);
t288 = r_i_i_C(2) * t240;
t241 = cos(t245);
t289 = r_i_i_C(1) * t241;
t260 = t234 - t288 + t289;
t271 = t286 * t251;
t297 = t260 * t254 + t271;
t296 = r_i_i_C(1) * t240 + r_i_i_C(2) * t241;
t295 = t244 * (t247 * t251 + qJD(1));
t250 = sin(qJ(4));
t243 = sin(t249);
t291 = pkin(5) * t243;
t231 = -t247 * t291 - t250 * t285;
t242 = qJD(6) + t247;
t258 = t296 * t242 - t231;
t252 = sin(qJ(1));
t281 = qJD(1) * t251;
t267 = t242 + t281;
t255 = cos(qJ(1));
t276 = qJD(3) * t255;
t272 = t254 * t276;
t293 = t267 * t252 - t272;
t277 = qJD(3) * t254;
t273 = t252 * t277;
t292 = t267 * t255 + t273;
t268 = -t242 * t251 - qJD(1);
t261 = t268 * t255;
t227 = t293 * t240 + t241 * t261;
t228 = t240 * t261 - t293 * t241;
t283 = -t227 * r_i_i_C(1) + t228 * r_i_i_C(2);
t262 = t268 * t252;
t229 = -t292 * t240 + t241 * t262;
t230 = t240 * t262 + t292 * t241;
t282 = t229 * r_i_i_C(1) - t230 * r_i_i_C(2);
t280 = qJD(1) * t252;
t279 = qJD(1) * t255;
t278 = qJD(3) * t251;
t275 = t242 * t289;
t274 = t254 * t242 * t288 + t296 * t278;
t265 = -t247 - t281;
t237 = t250 * pkin(4) + t291;
t263 = -t237 * t281 + t231;
t259 = -t254 * t275 + t274;
t257 = qJD(1) * t238 + t232 * t251 + t237 * t277;
t256 = t251 * t231 + qJD(2) + (t234 * t254 + t271) * qJD(3) + (-pkin(1) - pkin(7) - t237) * qJD(1);
t1 = [t228 * r_i_i_C(1) + t227 * r_i_i_C(2) + t299 * t252 + t256 * t255, t279, t297 * t279 + (-t258 * t254 + (-t260 * t251 + t270) * qJD(3)) * t252, -t257 * t252 + t263 * t255 + t282 (-t252 * t295 + (t265 * t255 - t273) * t243) * pkin(5) + t282, t282; t230 * r_i_i_C(1) + t229 * r_i_i_C(2) + t256 * t252 - t299 * t255, t280 (t260 * t276 + t286 * t280) * t251 + (t260 * t280 + (-t286 * qJD(3) + t258) * t255) * t254, t263 * t252 + t257 * t255 + t283 (t255 * t295 + (t265 * t252 + t272) * t243) * pkin(5) + t283, t283; 0, 0, -t297 * qJD(3) + t258 * t251, t237 * t278 + (-t232 - t275) * t254 + t274 (t243 * t278 - t254 * t284) * pkin(5) + t259, t259;];
JaD_transl  = t1;
