% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:42:35
% EndTime: 2019-02-26 22:42:36
% DurationCPUTime: 0.30s
% Computational Cost: add. (578->73), mult. (572->105), div. (0->0), fcn. (452->10), ass. (0->63)
t249 = qJ(3) + qJ(4);
t243 = sin(t249);
t250 = sin(qJ(3));
t284 = pkin(3) * qJD(3);
t247 = qJD(3) + qJD(4);
t290 = pkin(4) * t247;
t233 = -t243 * t290 - t250 * t284;
t244 = cos(t249);
t253 = cos(qJ(3));
t238 = t253 * pkin(3) + pkin(4) * t244;
t236 = pkin(2) + t238;
t251 = sin(qJ(2));
t254 = cos(qJ(2));
t285 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8);
t269 = t285 * t254;
t296 = (-t236 * t251 + t269) * qJD(2) + t254 * t233;
t255 = cos(qJ(1));
t242 = qJD(5) + t247;
t268 = t242 * t254 - qJD(1);
t294 = t255 * t268;
t245 = qJ(5) + t249;
t240 = sin(t245);
t241 = cos(t245);
t287 = r_i_i_C(2) * t241;
t265 = r_i_i_C(1) * t240 + t287;
t293 = -t242 * t265 + t233;
t277 = qJD(1) * t254;
t267 = -t242 + t277;
t252 = sin(qJ(1));
t275 = qJD(2) * t251;
t271 = t252 * t275;
t292 = t255 * t267 - t271;
t291 = pkin(4) * t243;
t289 = r_i_i_C(1) * t241;
t288 = r_i_i_C(2) * t240;
t237 = t250 * pkin(3) + t291;
t286 = pkin(7) + t237;
t282 = t242 * t251;
t270 = t255 * t275;
t257 = t252 * t267 + t270;
t229 = t240 * t257 - t241 * t294;
t230 = t240 * t294 + t241 * t257;
t280 = t229 * r_i_i_C(1) + t230 * r_i_i_C(2);
t262 = t268 * t252;
t231 = t292 * t240 + t241 * t262;
t232 = t240 * t262 - t292 * t241;
t279 = -t231 * r_i_i_C(1) + t232 * r_i_i_C(2);
t278 = qJD(1) * t252;
t276 = qJD(1) * t255;
t274 = qJD(2) * t254;
t273 = t244 * t290;
t272 = t242 * t289;
t266 = -t247 + t277;
t264 = t237 * t277 + t233;
t263 = t244 * (-t247 * t254 + qJD(1));
t261 = t236 - t288 + t289;
t260 = -t236 * t254 - t251 * t285 - pkin(1);
t259 = qJD(2) * t261;
t234 = t253 * t284 + t273;
t258 = qJD(1) * t238 - t234 * t254 + t237 * t275;
t256 = -qJD(2) * t285 - t293;
t235 = t282 * t288;
t1 = [t232 * r_i_i_C(1) + t231 * r_i_i_C(2) + t255 * t234 - t296 * t252 + (-t252 * t286 + t255 * t260) * qJD(1) (-t255 * t259 - t278 * t285) * t254 + (t255 * t256 + t261 * t278) * t251, t252 * t264 + t255 * t258 + t280 (t255 * t263 + (t252 * t266 + t270) * t243) * pkin(4) + t280, t280, 0; -t230 * r_i_i_C(1) + t229 * r_i_i_C(2) + t252 * t234 + t296 * t255 + (t252 * t260 + t255 * t286) * qJD(1) (-t252 * t259 + t276 * t285) * t254 + (t252 * t256 - t261 * t276) * t251, t252 * t258 - t255 * t264 + t279 (t252 * t263 + (-t255 * t266 + t271) * t243) * pkin(4) + t279, t279, 0; 0, t293 * t254 + (-t251 * t261 + t269) * qJD(2), t235 + (-t234 - t272) * t251 + (-t237 - t265) * t274, t235 + (-t272 - t273) * t251 + (-t265 - t291) * t274, -t274 * t287 + t235 + (-t240 * t274 - t241 * t282) * r_i_i_C(1), 0;];
JaD_transl  = t1;
