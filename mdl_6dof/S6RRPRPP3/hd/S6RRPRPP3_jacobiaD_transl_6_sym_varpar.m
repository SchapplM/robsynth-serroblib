% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPP3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:35:57
% EndTime: 2019-02-26 21:35:57
% DurationCPUTime: 0.35s
% Computational Cost: add. (449->66), mult. (717->99), div. (0->0), fcn. (613->8), ass. (0->48)
t251 = cos(pkin(9)) * pkin(3) + pkin(2);
t257 = sin(qJ(2));
t259 = cos(qJ(2));
t278 = pkin(5) + r_i_i_C(1) + pkin(8) + qJ(3);
t296 = t278 * t259;
t301 = (-t251 * t257 + t296) * qJD(2) + t257 * qJD(3);
t254 = pkin(9) + qJ(4);
t252 = sin(t254);
t253 = cos(t254);
t279 = pkin(4) + r_i_i_C(3) + qJ(6);
t293 = r_i_i_C(2) + qJ(5);
t294 = t279 * t252 - t253 * t293;
t300 = -t278 * qJD(2) + t294 * qJD(4) - qJD(5) * t252 - qJD(6) * t253;
t265 = -t252 * t293 - t279 * t253;
t263 = -t251 + t265;
t298 = t257 * t263 + t296;
t258 = sin(qJ(1));
t291 = t258 * t259;
t260 = cos(qJ(1));
t290 = t260 * t252;
t289 = t260 * t253;
t288 = qJD(1) * t258;
t287 = qJD(1) * t260;
t286 = qJD(2) * t257;
t285 = qJD(2) * t259;
t284 = qJD(2) * t260;
t283 = qJD(4) * t257;
t282 = qJD(4) * t258;
t281 = qJD(4) * t260;
t277 = pkin(3) * sin(pkin(9)) + pkin(7);
t276 = t258 * t286;
t275 = t257 * t284;
t274 = t252 * t282;
t273 = t252 * t281;
t272 = t253 * t281;
t268 = t258 * t252 + t259 * t289;
t237 = t252 * t291 + t289;
t267 = t252 * t287 + t253 * t282;
t266 = -t251 * t259 - t257 * t278 - pkin(1);
t262 = qJD(2) * t263 + qJD(3);
t261 = t300 * t257 + t262 * t259;
t239 = -t258 * t253 + t259 * t290;
t238 = t253 * t291 - t290;
t236 = qJD(1) * t268 - t253 * t276 - t259 * t274 - t272;
t235 = -t252 * t276 - t253 * t288 + t259 * t267 - t273;
t234 = t259 * t273 + (t259 * t288 + t275) * t253 - t267;
t233 = qJD(1) * t237 + t252 * t275 - t259 * t272 - t274;
t1 = [-t237 * qJD(5) - t238 * qJD(6) - t293 * t235 - t279 * t236 - t301 * t258 + (-t258 * t277 + t260 * t266) * qJD(1), t261 * t260 - t298 * t288, -t257 * t288 + t259 * t284, qJD(5) * t268 - t239 * qJD(6) + t279 * t233 - t234 * t293, -t233, -t234; t239 * qJD(5) + t268 * qJD(6) - t293 * t233 - t279 * t234 + t301 * t260 + (t258 * t266 + t260 * t277) * qJD(1), t261 * t258 + t298 * t287, t257 * t287 + t258 * t285, t238 * qJD(5) - t237 * qJD(6) - t279 * t235 + t236 * t293, t235, t236; 0, t262 * t257 - t300 * t259, t286, -t294 * t285 + (qJD(4) * t265 + qJD(5) * t253 - qJD(6) * t252) * t257, t252 * t285 + t253 * t283, -t252 * t283 + t253 * t285;];
JaD_transl  = t1;
