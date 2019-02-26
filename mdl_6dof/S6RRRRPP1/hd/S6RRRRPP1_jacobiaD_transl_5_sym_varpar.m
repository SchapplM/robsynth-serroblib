% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:25:04
% EndTime: 2019-02-26 22:25:04
% DurationCPUTime: 0.35s
% Computational Cost: add. (443->75), mult. (515->107), div. (0->0), fcn. (394->10), ass. (0->63)
t264 = qJ(4) + pkin(10);
t259 = sin(t264);
t260 = cos(t264);
t330 = r_i_i_C(1) * t259 + r_i_i_C(2) * t260;
t263 = qJD(2) + qJD(3);
t317 = r_i_i_C(2) * t259;
t329 = t263 * t317 + qJD(5);
t270 = cos(qJ(4));
t314 = t270 * pkin(4);
t257 = pkin(3) + t314;
t318 = r_i_i_C(1) * t260;
t328 = t257 + t318;
t265 = qJ(2) + qJ(3);
t261 = sin(t265);
t262 = cos(t265);
t268 = sin(qJ(2));
t312 = pkin(2) * qJD(2);
t300 = t268 * t312;
t267 = sin(qJ(4));
t311 = pkin(4) * qJD(4);
t301 = t267 * t311;
t266 = -qJ(5) - pkin(9);
t313 = r_i_i_C(3) - t266;
t320 = pkin(4) * t267;
t327 = (t313 * t263 - t301) * t262 + (pkin(8) + pkin(7) + t320) * qJD(1) - (t257 * t263 - qJD(5)) * t261 - t300;
t304 = qJD(4) * t261;
t326 = t261 * t301 + t329 * t262 + t330 * t304;
t272 = cos(qJ(1));
t288 = qJD(4) * t262 - qJD(1);
t325 = t272 * t288;
t287 = qJD(1) * t262 - qJD(4);
t269 = sin(qJ(1));
t308 = t263 * t269;
t299 = t261 * t308;
t323 = t287 * t272 - t299;
t321 = pkin(2) * t268;
t315 = r_i_i_C(3) * t262;
t310 = t262 * t263;
t309 = t262 * t266;
t307 = t263 * t272;
t306 = qJD(1) * t269;
t305 = qJD(1) * t272;
t298 = t261 * t307;
t297 = t261 * t306;
t296 = t261 * t305;
t286 = t328 * t263;
t285 = t288 * t270;
t284 = t288 * t269;
t283 = -t320 - t330;
t282 = -t261 * t328 - t309;
t281 = t266 * t299 + t326 * t269 + t296 * t317 + t305 * t315;
t280 = t266 * t298 + t326 * t272 + t328 * t297 + t306 * t309;
t279 = (-r_i_i_C(3) * t261 - t262 * t328) * t263;
t277 = t287 * t269 + t298;
t271 = cos(qJ(2));
t276 = t270 * t311 + (-t271 * pkin(2) - t257 * t262 - t313 * t261 - pkin(1)) * qJD(1);
t275 = -t271 * t312 + t279;
t274 = r_i_i_C(3) * t310 + (t283 * qJD(4) - t263 * t266) * t262 + (-t286 + t329) * t261;
t236 = t259 * t284 - t260 * t323;
t235 = t323 * t259 + t260 * t284;
t234 = t259 * t325 + t277 * t260;
t233 = t277 * t259 - t260 * t325;
t1 = [t236 * r_i_i_C(1) + t235 * r_i_i_C(2) - t327 * t269 + t276 * t272 (-t261 * t317 - t315 + t321) * t306 + t275 * t272 + t280 (-r_i_i_C(3) * t307 - t306 * t317) * t261 + (-r_i_i_C(3) * t306 - t272 * t286) * t262 + t280, t233 * r_i_i_C(1) + t234 * r_i_i_C(2) + (t277 * t267 - t272 * t285) * pkin(4), t262 * t307 - t297, 0; -t234 * r_i_i_C(1) + t233 * r_i_i_C(2) + t276 * t269 + t327 * t272, t275 * t269 + (t282 - t321) * t305 + t281, t269 * t279 + t282 * t305 + t281, -t235 * r_i_i_C(1) + t236 * r_i_i_C(2) + (-t267 * t323 - t269 * t285) * pkin(4), t262 * t308 + t296, 0; 0, t274 - t300, t274, t283 * t310 + (-t314 + t317 - t318) * t304, t263 * t261, 0;];
JaD_transl  = t1;
