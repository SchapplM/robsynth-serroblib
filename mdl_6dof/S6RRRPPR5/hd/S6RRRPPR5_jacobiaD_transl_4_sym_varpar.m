% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR5_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR5_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobiaD_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:46
% EndTime: 2019-02-26 22:05:47
% DurationCPUTime: 0.26s
% Computational Cost: add. (252->77), mult. (593->125), div. (0->0), fcn. (554->10), ass. (0->52)
t253 = qJ(3) + pkin(11);
t251 = sin(t253);
t252 = cos(t253);
t257 = sin(qJ(3));
t291 = pkin(3) * t257;
t264 = r_i_i_C(1) * t251 + r_i_i_C(2) * t252 + t291;
t263 = qJD(3) * t264;
t290 = t251 * r_i_i_C(2);
t260 = cos(qJ(3));
t289 = t260 * pkin(3);
t288 = r_i_i_C(3) + qJ(4) + pkin(9);
t254 = sin(pkin(6));
t258 = sin(qJ(2));
t287 = t254 * t258;
t259 = sin(qJ(1));
t286 = t254 * t259;
t285 = t254 * t260;
t262 = cos(qJ(1));
t284 = t254 * t262;
t283 = t259 * t258;
t261 = cos(qJ(2));
t282 = t259 * t261;
t281 = t262 * t258;
t280 = t262 * t261;
t279 = qJD(1) * t259;
t278 = qJD(1) * t262;
t277 = qJD(2) * t258;
t276 = qJD(2) * t261;
t255 = cos(pkin(6));
t241 = t255 * t281 + t282;
t275 = qJD(3) * t241;
t274 = qJD(3) * t262;
t273 = t255 * t283;
t272 = t255 * t280;
t271 = t254 * t279;
t270 = t254 * t278;
t269 = t254 * t274;
t268 = qJD(2) * t255 + qJD(1);
t250 = pkin(2) + t289;
t267 = t252 * r_i_i_C(1) + t250 - t290;
t266 = t255 * t282 + t281;
t265 = t271 - t275;
t244 = t252 * t269;
t243 = -t273 + t280;
t240 = -t272 + t283;
t239 = -qJD(1) * t273 - t259 * t277 + t268 * t280;
t238 = t266 * qJD(1) + t241 * qJD(2);
t237 = t241 * qJD(1) + t266 * qJD(2);
t236 = -qJD(1) * t272 - t262 * t276 + t268 * t283;
t235 = t251 * t270 - t237 * t252 + (-t243 * t251 + t252 * t286) * qJD(3);
t234 = t252 * t270 + t237 * t251 + (-t243 * t252 - t251 * t286) * qJD(3);
t1 = [(-t239 * t252 + t251 * t275 + t244) * r_i_i_C(1) + (t239 * t251 + t252 * t275) * r_i_i_C(2) - t239 * t250 + t275 * t291 - t240 * qJD(4) - pkin(1) * t278 - t288 * t238 + ((t289 - t290) * t274 + (-pkin(8) - t264) * t279) * t254, t243 * qJD(4) + t267 * t236 - t288 * t237 + t263 * t266, t234 * r_i_i_C(1) - t235 * r_i_i_C(2) + (t260 * t270 + t237 * t257 + (-t243 * t260 - t257 * t286) * qJD(3)) * pkin(3), -t236, 0, 0; t235 * r_i_i_C(1) + t234 * r_i_i_C(2) + t266 * qJD(4) - t237 * t250 - t288 * t236 + (-pkin(1) * t259 + pkin(8) * t284) * qJD(1) + (t257 * t270 + (-t243 * t257 + t259 * t285) * qJD(3)) * pkin(3), t241 * qJD(4) - t267 * t238 + t288 * t239 + t240 * t263, t244 * r_i_i_C(2) + (t265 * r_i_i_C(1) - t239 * r_i_i_C(2)) * t252 + ((-t239 + t269) * r_i_i_C(1) - t265 * r_i_i_C(2)) * t251 + (t260 * t271 - t239 * t257 + (-t241 * t260 + t257 * t284) * qJD(3)) * pkin(3), t238, 0, 0; 0 (qJD(4) * t258 - t261 * t263 + (-t267 * t258 + t288 * t261) * qJD(2)) * t254, -t264 * t254 * t276 + ((-t251 * t255 - t252 * t287) * r_i_i_C(1) + (t251 * t287 - t252 * t255) * r_i_i_C(2) + (-t255 * t257 - t258 * t285) * pkin(3)) * qJD(3), t254 * t277, 0, 0;];
JaD_transl  = t1;
