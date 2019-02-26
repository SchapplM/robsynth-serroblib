% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:56:18
% EndTime: 2019-02-26 20:56:19
% DurationCPUTime: 0.40s
% Computational Cost: add. (534->71), mult. (614->104), div. (0->0), fcn. (508->10), ass. (0->55)
t270 = cos(qJ(4));
t309 = pkin(4) * t270;
t260 = pkin(3) + t309;
t268 = sin(qJ(4));
t269 = sin(qJ(3));
t271 = cos(qJ(3));
t265 = qJ(4) + pkin(10);
t261 = sin(t265);
t308 = r_i_i_C(2) + qJ(5) + pkin(8);
t280 = t308 * qJD(3) + qJD(6) * t261;
t306 = pkin(4) * qJD(4);
t320 = (-t268 * t306 + t280) * t271 - (qJD(3) * t260 - qJD(5)) * t269;
t263 = cos(t265);
t307 = r_i_i_C(3) + qJ(6);
t310 = pkin(4) * t268;
t311 = pkin(5) + r_i_i_C(1);
t312 = t311 * t261 - t307 * t263 + t310;
t318 = qJD(4) * t312 - t280;
t278 = -t307 * t261 - t311 * t263;
t276 = -t260 + t278;
t317 = t276 * t269 + t308 * t271;
t299 = qJD(1) * t271;
t316 = t268 * (-qJD(4) + t299);
t266 = qJ(1) + pkin(9);
t262 = sin(t266);
t305 = t262 * t261;
t304 = t262 * t271;
t264 = cos(t266);
t303 = t264 * t263;
t302 = qJD(1) * t262;
t301 = qJD(1) * t264;
t300 = qJD(1) * t269;
t298 = qJD(3) * t269;
t297 = qJD(3) * t271;
t296 = qJD(4) * t263;
t295 = qJD(4) * t264;
t293 = qJD(6) * t263;
t291 = pkin(7) + t310;
t290 = t264 * t298;
t289 = t262 * t298;
t288 = qJD(4) * t305;
t287 = t261 * t295;
t286 = t263 * t295;
t283 = t270 * t306 - t293;
t282 = t271 * t303 + t305;
t279 = -t260 * t271 - t308 * t269 - pkin(2);
t277 = t261 * t301 + t262 * t296;
t275 = t268 * t298 + (-qJD(4) * t271 + qJD(1)) * t270;
t273 = t276 * qJD(3) + qJD(5);
t272 = t318 * t269 + t273 * t271;
t249 = t282 * qJD(1) - t263 * t289 - t271 * t288 - t286;
t248 = -t261 * t289 - t263 * t302 + t277 * t271 - t287;
t247 = t271 * t287 + (t262 * t299 + t290) * t263 - t277;
t246 = t261 * t290 - t271 * t286 - t288 + (t261 * t304 + t303) * qJD(1);
t1 = [t283 * t264 - t311 * t249 - t307 * t248 - t320 * t262 + (-cos(qJ(1)) * pkin(1) - t291 * t262 + t279 * t264) * qJD(1), 0, t272 * t264 - t317 * t302, t282 * qJD(6) - t307 * t247 + t311 * t246 + (t262 * t316 + t275 * t264) * pkin(4), -t262 * t300 + t264 * t297, -t246; t283 * t262 - t311 * t247 - t307 * t246 + t320 * t264 + (-sin(qJ(1)) * pkin(1) + t291 * t264 + t279 * t262) * qJD(1), 0, t272 * t262 + t317 * t301 -(t264 * t261 - t263 * t304) * qJD(6) + t307 * t249 - t311 * t248 + (t275 * t262 - t264 * t316) * pkin(4), t262 * t297 + t264 * t300, t248; 0, 0, t273 * t269 - t318 * t271, -t312 * t297 + (t293 + (t278 - t309) * qJD(4)) * t269, t298, t261 * t297 + t269 * t296;];
JaD_transl  = t1;
