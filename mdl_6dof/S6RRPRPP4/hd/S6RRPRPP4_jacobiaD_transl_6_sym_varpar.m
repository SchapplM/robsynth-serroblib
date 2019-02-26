% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPP4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:36:37
% EndTime: 2019-02-26 21:36:37
% DurationCPUTime: 0.26s
% Computational Cost: add. (372->67), mult. (680->96), div. (0->0), fcn. (558->8), ass. (0->51)
t249 = sin(qJ(2));
t251 = cos(qJ(4));
t252 = cos(qJ(2));
t246 = qJ(4) + pkin(9);
t245 = cos(t246);
t275 = pkin(2) + r_i_i_C(2) + qJ(5) + pkin(8);
t262 = -t275 * qJD(2) - t245 * qJD(6) + qJD(3);
t248 = sin(qJ(4));
t289 = pkin(4) * t248;
t270 = qJ(3) + t289;
t285 = pkin(4) * qJD(4);
t288 = t251 * pkin(4);
t297 = (t251 * t285 + t262) * t249 + (t270 * qJD(2) + qJD(5)) * t252 + (pkin(7) + pkin(3) + t288) * qJD(1);
t296 = t275 * t249;
t250 = sin(qJ(1));
t278 = qJD(2) * t252;
t272 = t250 * t278;
t253 = cos(qJ(1));
t280 = qJD(1) * t253;
t294 = t249 * t280 + t272;
t244 = sin(t246);
t286 = r_i_i_C(3) + qJ(6);
t290 = -r_i_i_C(1) - pkin(5);
t291 = t290 * t244 + t286 * t245;
t284 = t250 * t245;
t283 = t253 * t244;
t282 = t253 * t245;
t281 = qJD(1) * t250;
t279 = qJD(2) * t249;
t277 = qJD(2) * t253;
t276 = t244 * qJD(6);
t271 = t252 * t277;
t269 = qJD(4) * t249 + qJD(1);
t268 = qJD(1) * t249 + qJD(4);
t266 = t269 * t248;
t265 = t269 * t250;
t264 = t268 * t253;
t263 = t249 * t283 + t284;
t261 = -t268 * t250 + t271;
t260 = t289 - t291;
t259 = t286 * t244 - t290 * t245 + t288;
t258 = -t270 + t291;
t257 = qJ(3) + t260;
t256 = -t248 * t285 + t276 + (-t270 * t249 - t275 * t252 - pkin(1)) * qJD(1);
t255 = t259 * qJD(4) + t262;
t254 = (t258 * qJD(2) - qJD(5)) * t249 + t255 * t252;
t238 = t261 * t244 + t269 * t282;
t237 = -t245 * t271 + t263 * qJD(4) + (t249 * t284 + t283) * qJD(1);
t236 = t245 * t265 + (t264 + t272) * t244;
t235 = -qJD(4) * t282 + t244 * t265 - t294 * t245;
t1 = [-t286 * t235 + t290 * t236 - t297 * t250 + t256 * t253 (t258 * t252 + t296) * t281 + t254 * t253, -t249 * t281 + t271, t263 * qJD(6) + t286 * t238 + t290 * t237 + (t261 * t251 - t253 * t266) * pkin(4), -t249 * t277 - t252 * t281, t237; t286 * t237 - t290 * t238 + t256 * t250 + t297 * t253 (t257 * t252 - t296) * t280 + t254 * t250, t294 -(-t250 * t249 * t244 + t282) * qJD(6) + t286 * t236 + t290 * t235 + (t251 * t264 + (t251 * t278 - t266) * t250) * pkin(4), -t250 * t279 + t252 * t280, t235; 0 (t257 * qJD(2) + qJD(5)) * t252 + t255 * t249, t279, t259 * t279 + (t260 * qJD(4) - t276) * t252, t278, -t252 * qJD(4) * t244 - t245 * t279;];
JaD_transl  = t1;
