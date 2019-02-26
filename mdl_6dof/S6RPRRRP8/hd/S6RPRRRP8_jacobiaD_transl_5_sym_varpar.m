% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP8_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP8_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:11:55
% EndTime: 2019-02-26 21:11:55
% DurationCPUTime: 0.26s
% Computational Cost: add. (274->62), mult. (406->96), div. (0->0), fcn. (305->8), ass. (0->56)
t302 = pkin(9) + r_i_i_C(3);
t258 = cos(qJ(5));
t299 = r_i_i_C(1) * t258;
t308 = pkin(4) + t299;
t256 = sin(qJ(3));
t254 = qJ(3) + qJ(4);
t252 = cos(t254);
t285 = t302 * t252;
t251 = sin(t254);
t301 = pkin(4) * t251;
t306 = -pkin(3) * t256 - qJ(2) + t285 - t301;
t260 = cos(qJ(1));
t271 = qJD(1) * t251 + qJD(5);
t305 = t271 * t260;
t304 = t302 * t251;
t257 = sin(qJ(1));
t253 = qJD(3) + qJD(4);
t293 = t253 * t260;
t303 = -t252 * t293 + t271 * t257;
t300 = pkin(4) * t252;
t255 = sin(qJ(5));
t298 = r_i_i_C(2) * t255;
t297 = -pkin(1) - pkin(8) - pkin(7);
t296 = pkin(3) * qJD(3);
t295 = t251 * t255;
t294 = t253 * t258;
t292 = qJD(1) * t257;
t291 = qJD(1) * t260;
t290 = qJD(5) * t255;
t289 = qJD(5) * t258;
t288 = t252 * t298;
t287 = t256 * t296;
t259 = cos(qJ(3));
t286 = t259 * t296;
t284 = t253 * t295;
t282 = t252 * t253 * t257;
t278 = t252 * t291;
t277 = t252 * t290;
t276 = t252 * t289;
t275 = r_i_i_C(2) * t284;
t274 = qJD(1) * t252 * t299;
t273 = r_i_i_C(2) * t276;
t272 = -qJD(5) * t251 - qJD(1);
t270 = t272 * t260;
t268 = pkin(4) * t278 + t257 * t275 + t260 * t274 + t302 * (t251 * t291 + t282);
t267 = qJD(1) * (pkin(3) * t259 - t288);
t266 = t251 * t294 + t277;
t265 = t257 * t274 + t292 * t300 + (r_i_i_C(1) * t277 + t273) * t260 + (t302 * t292 + t308 * t293) * t251;
t264 = t286 + qJD(2) + (t300 + t304) * t253;
t263 = (-t308 * t252 + t288 - t304) * t253 + (r_i_i_C(1) * t290 + r_i_i_C(2) * t289) * t251;
t262 = -t266 * r_i_i_C(1) - t253 * t301 - t273;
t232 = t258 * t305 + (t252 * t294 + t272 * t255) * t257;
t231 = t272 * t258 * t257 + (-t282 - t305) * t255;
t230 = t255 * t270 - t303 * t258;
t229 = t303 * t255 + t258 * t270;
t1 = [t230 * r_i_i_C(1) + t229 * r_i_i_C(2) + t264 * t260 + (t306 * t257 + t297 * t260) * qJD(1), t291, t260 * t267 + (t262 - t287) * t257 + t268, t262 * t257 - t278 * t298 + t268, t231 * r_i_i_C(1) - t232 * r_i_i_C(2), 0; t232 * r_i_i_C(1) + t231 * r_i_i_C(2) + t264 * t257 + (t297 * t257 - t306 * t260) * qJD(1), t292, t257 * t267 + (t287 + (-r_i_i_C(2) * t295 - t285) * t253) * t260 + t265, -t260 * t275 + (-t292 * t298 - t302 * t293) * t252 + t265, -t229 * r_i_i_C(1) + t230 * r_i_i_C(2), 0; 0, 0, t263 - t286, t263, t266 * r_i_i_C(2) + (-t276 + t284) * r_i_i_C(1), 0;];
JaD_transl  = t1;
