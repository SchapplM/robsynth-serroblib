% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:47:46
% EndTime: 2019-02-26 21:47:47
% DurationCPUTime: 0.32s
% Computational Cost: add. (374->56), mult. (446->81), div. (0->0), fcn. (351->10), ass. (0->52)
t255 = cos(qJ(4));
t242 = t255 * pkin(4) + pkin(3);
t249 = qJ(2) + pkin(10);
t244 = sin(t249);
t245 = cos(t249);
t289 = r_i_i_C(3) + pkin(9) + pkin(8);
t264 = t289 * t245 - sin(qJ(2)) * pkin(2);
t270 = qJD(4) * t245 - qJD(1);
t252 = sin(qJ(4));
t293 = pkin(4) * t252;
t302 = (-t242 * t244 + t264) * qJD(2) - t270 * t293 - qJD(1) * (-qJ(3) - pkin(7));
t250 = qJ(4) + qJ(5);
t247 = cos(t250);
t246 = sin(t250);
t292 = r_i_i_C(2) * t246;
t265 = r_i_i_C(1) * t247 + t242 - t292;
t261 = -t265 * t244 + t264;
t257 = cos(qJ(1));
t248 = qJD(4) + qJD(5);
t272 = t245 * t248 - qJD(1);
t300 = t257 * t272;
t298 = -t289 * t244 - cos(qJ(2)) * pkin(2);
t291 = r_i_i_C(2) * t247;
t268 = r_i_i_C(1) * t246 + t291;
t288 = pkin(4) * qJD(4);
t297 = t268 * t248 + t252 * t288;
t283 = qJD(1) * t245;
t271 = -t248 + t283;
t254 = sin(qJ(1));
t279 = qJD(2) * t244;
t275 = t254 * t279;
t296 = t271 * t257 - t275;
t286 = t247 * t248;
t274 = t257 * t279;
t262 = t271 * t254 + t274;
t237 = t262 * t246 - t247 * t300;
t238 = t246 * t300 + t262 * t247;
t285 = t237 * r_i_i_C(1) + t238 * r_i_i_C(2);
t267 = t272 * t254;
t239 = t296 * t246 + t247 * t267;
t240 = t246 * t267 - t296 * t247;
t284 = -t239 * r_i_i_C(1) + t240 * r_i_i_C(2);
t281 = qJD(1) * t254;
t280 = qJD(1) * t257;
t278 = qJD(2) * t245;
t276 = t255 * t288;
t269 = -qJD(4) + t283;
t266 = t270 * t255;
t260 = t276 + qJD(3) + (-t242 * t245 - pkin(1) + t298) * qJD(1);
t259 = t297 * t244 + (-t265 * t245 + t298) * qJD(2);
t241 = t244 * t248 * t292;
t1 = [t240 * r_i_i_C(1) + t239 * r_i_i_C(2) - t302 * t254 + t260 * t257, t259 * t257 - t261 * t281, t280 (-t257 * t266 + (t269 * t254 + t274) * t252) * pkin(4) + t285, t285, 0; -t238 * r_i_i_C(1) + t237 * r_i_i_C(2) + t260 * t254 + t302 * t257, t259 * t254 + t261 * t280, t281 (-t254 * t266 + (-t269 * t257 + t275) * t252) * pkin(4) + t284, t284, 0; 0, t261 * qJD(2) - t297 * t245, 0, t241 + (-r_i_i_C(1) * t286 - t276) * t244 + (-t268 - t293) * t278, -t278 * t291 + t241 + (-t244 * t286 - t246 * t278) * r_i_i_C(1), 0;];
JaD_transl  = t1;
