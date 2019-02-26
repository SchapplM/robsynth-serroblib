% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:17:37
% EndTime: 2019-02-26 21:17:37
% DurationCPUTime: 0.28s
% Computational Cost: add. (367->60), mult. (426->91), div. (0->0), fcn. (338->9), ass. (0->53)
t251 = cos(qJ(4));
t240 = pkin(4) * t251 + pkin(3);
t245 = pkin(11) + qJ(3);
t241 = sin(t245);
t242 = cos(t245);
t283 = r_i_i_C(3) + pkin(9) + pkin(8);
t268 = t283 * t242;
t249 = sin(qJ(4));
t282 = pkin(4) * qJD(4);
t272 = t249 * t282;
t291 = (-t240 * t241 + t268) * qJD(3) - t242 * t272;
t252 = cos(qJ(1));
t246 = qJD(4) + qJD(5);
t266 = t242 * t246 - qJD(1);
t289 = t252 * t266;
t247 = qJ(4) + qJ(5);
t243 = sin(t247);
t244 = cos(t247);
t284 = r_i_i_C(2) * t244;
t261 = r_i_i_C(1) * t243 + t284;
t288 = -t261 * t246 - t272;
t277 = qJD(1) * t242;
t265 = -t246 + t277;
t250 = sin(qJ(1));
t274 = qJD(3) * t241;
t270 = t250 * t274;
t287 = t265 * t252 - t270;
t286 = pkin(4) * t249;
t285 = r_i_i_C(2) * t243;
t280 = t244 * t246;
t269 = t252 * t274;
t255 = t265 * t250 + t269;
t234 = t255 * t243 - t244 * t289;
t235 = t243 * t289 + t255 * t244;
t279 = t234 * r_i_i_C(1) + t235 * r_i_i_C(2);
t260 = t266 * t250;
t236 = t287 * t243 + t244 * t260;
t237 = t243 * t260 - t287 * t244;
t278 = -t236 * r_i_i_C(1) + t237 * r_i_i_C(2);
t276 = qJD(1) * t250;
t275 = qJD(1) * t252;
t273 = qJD(3) * t242;
t271 = t251 * t282;
t267 = pkin(7) + qJ(2) + t286;
t263 = -qJD(4) + t277;
t262 = qJD(2) + t271;
t259 = (-qJD(4) * t242 + qJD(1)) * t251;
t258 = r_i_i_C(1) * t244 + t240 - t285;
t257 = qJD(3) * t258;
t256 = -t240 * t242 - t283 * t241 - cos(pkin(11)) * pkin(2) - pkin(1);
t254 = -t283 * qJD(3) - t288;
t238 = t241 * t246 * t285;
t1 = [t237 * r_i_i_C(1) + t236 * r_i_i_C(2) + t262 * t252 - t291 * t250 + (-t267 * t250 + t256 * t252) * qJD(1), t275 (-t252 * t257 - t283 * t276) * t242 + (t254 * t252 + t258 * t276) * t241 (t252 * t259 + (t263 * t250 + t269) * t249) * pkin(4) + t279, t279, 0; -t235 * r_i_i_C(1) + t234 * r_i_i_C(2) + t262 * t250 + t291 * t252 + (t256 * t250 + t267 * t252) * qJD(1), t276 (-t250 * t257 + t283 * t275) * t242 + (t254 * t250 - t258 * t275) * t241 (t250 * t259 + (-t263 * t252 + t270) * t249) * pkin(4) + t278, t278, 0; 0, 0, t288 * t242 + (-t258 * t241 + t268) * qJD(3), t238 + (-r_i_i_C(1) * t280 - t271) * t241 + (-t261 - t286) * t273, -t273 * t284 + t238 + (-t241 * t280 - t243 * t273) * r_i_i_C(1), 0;];
JaD_transl  = t1;
