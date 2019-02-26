% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function JaD_transl = S6RRPRPP3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:35:57
% EndTime: 2019-02-26 21:35:57
% DurationCPUTime: 0.26s
% Computational Cost: add. (335->61), mult. (543->94), div. (0->0), fcn. (458->8), ass. (0->45)
t256 = sin(qJ(2));
t250 = cos(pkin(9)) * pkin(3) + pkin(2);
t253 = pkin(9) + qJ(4);
t251 = sin(t253);
t252 = cos(t253);
t290 = r_i_i_C(3) + qJ(5);
t292 = pkin(4) - r_i_i_C(2);
t293 = t290 * t251 + t292 * t252 + t250;
t258 = cos(qJ(2));
t291 = r_i_i_C(1) + pkin(8) + qJ(3);
t295 = t291 * t258;
t299 = t293 * t256 - t295;
t276 = t256 * qJD(3);
t278 = qJD(5) * t251;
t298 = (-t250 * t256 + t295) * qJD(2) + t258 * t278 + t276;
t296 = -t278 + (t292 * t251 - t290 * t252) * qJD(4);
t257 = sin(qJ(1));
t288 = t257 * t258;
t259 = cos(qJ(1));
t287 = t259 * t252;
t286 = qJD(1) * t257;
t285 = qJD(1) * t259;
t284 = qJD(2) * t256;
t283 = qJD(2) * t258;
t282 = qJD(2) * t259;
t281 = qJD(4) * t256;
t280 = qJD(4) * t257;
t279 = qJD(4) * t259;
t277 = t252 * qJD(5);
t275 = pkin(3) * sin(pkin(9)) + pkin(7);
t274 = t257 * t284;
t273 = t251 * t280;
t272 = t256 * t282;
t271 = t251 * t279;
t270 = t252 * t279;
t269 = t291 * t256;
t266 = t257 * t251 + t258 * t287;
t264 = -t250 * t258 - pkin(1) - t269;
t263 = t251 * t285 + t252 * t280;
t260 = qJD(3) * t258 + t296 * t256 + (-t258 * t293 - t269) * qJD(2);
t239 = t266 * qJD(1) - t252 * t274 - t258 * t273 - t270;
t238 = -t251 * t274 - t252 * t286 + t263 * t258 - t271;
t237 = t258 * t271 + (t258 * t286 + t272) * t252 - t263;
t236 = t251 * t272 - t258 * t270 - t273 + (t251 * t288 + t287) * qJD(1);
t1 = [-t259 * t277 - t292 * t239 - t290 * t238 - t298 * t257 + (-t275 * t257 + t264 * t259) * qJD(1), t260 * t259 + t299 * t286, -t256 * t286 + t258 * t282, t266 * qJD(5) + t292 * t236 - t290 * t237, -t236, 0; -t257 * t277 - t292 * t237 - t290 * t236 + t298 * t259 + (t264 * t257 + t275 * t259) * qJD(1), t260 * t257 - t285 * t299, t256 * t285 + t257 * t283 -(t259 * t251 - t252 * t288) * qJD(5) + t290 * t239 - t292 * t238, t238, 0; 0, -qJD(2) * t299 - t296 * t258 + t276, t284 (-t290 * t281 - t292 * t283) * t251 + (t290 * t283 + (-t292 * qJD(4) + qJD(5)) * t256) * t252, t251 * t283 + t252 * t281, 0;];
JaD_transl  = t1;
