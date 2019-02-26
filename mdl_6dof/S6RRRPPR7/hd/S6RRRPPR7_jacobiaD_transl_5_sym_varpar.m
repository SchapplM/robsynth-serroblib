% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:12
% EndTime: 2019-02-26 22:07:12
% DurationCPUTime: 0.31s
% Computational Cost: add. (271->64), mult. (846->97), div. (0->0), fcn. (752->8), ass. (0->50)
t262 = sin(qJ(2));
t265 = cos(qJ(2));
t290 = t262 * qJD(5);
t261 = sin(qJ(3));
t291 = t261 * qJD(4);
t288 = -r_i_i_C(3) - qJ(5) + pkin(8);
t304 = t288 * t265;
t307 = (-t262 * pkin(2) + t304) * qJD(2) + t265 * t291 - t290;
t264 = cos(qJ(3));
t260 = cos(pkin(10));
t259 = sin(pkin(10));
t281 = t259 * r_i_i_C(2) - pkin(3) - pkin(4);
t276 = t260 * r_i_i_C(1) - t281;
t282 = -t259 * r_i_i_C(1) - qJ(4);
t277 = t260 * r_i_i_C(2) - t282;
t303 = t276 * t261 - t277 * t264;
t306 = t303 * qJD(3) - t291;
t271 = -t277 * t261 - t276 * t264;
t269 = -pkin(2) + t271;
t268 = t269 * t262 + t304;
t263 = sin(qJ(1));
t301 = t263 * t261;
t300 = t263 * t265;
t266 = cos(qJ(1));
t299 = t266 * t264;
t298 = qJD(1) * t263;
t297 = qJD(1) * t266;
t296 = qJD(2) * t262;
t295 = qJD(2) * t265;
t294 = qJD(2) * t266;
t293 = qJD(3) * t264;
t292 = qJD(3) * t266;
t289 = t264 * qJD(4);
t287 = t263 * t296;
t286 = qJD(3) * t301;
t285 = t262 * t294;
t284 = t261 * t292;
t283 = t264 * t292;
t280 = t288 * t262;
t275 = t265 * t299 + t301;
t273 = t261 * t297 + t263 * t293;
t272 = -pkin(2) * t265 - pkin(1) - t280;
t267 = -t265 * qJD(5) + t306 * t262 + (t269 * t265 - t280) * qJD(2);
t248 = t275 * qJD(1) - t264 * t287 - t265 * t286 - t283;
t247 = -t261 * t287 - t264 * t298 + t273 * t265 - t284;
t246 = t265 * t284 + (t265 * t298 + t285) * t264 - t273;
t245 = t261 * t285 - t265 * t283 - t286 + (t261 * t300 + t299) * qJD(1);
t244 = t247 * t260;
t243 = t246 * t260;
t1 = [-t266 * t289 - t244 * r_i_i_C(2) + t282 * t247 - t276 * t248 - t307 * t263 + (-t263 * pkin(7) + t272 * t266) * qJD(1), t267 * t266 - t268 * t298, -t243 * r_i_i_C(2) + t275 * qJD(4) + t276 * t245 + t282 * t246, -t245, t262 * t298 - t265 * t294, 0; -t263 * t289 - t243 * r_i_i_C(1) - t277 * t245 + t281 * t246 + t307 * t266 + (t266 * pkin(7) + t272 * t263) * qJD(1), t267 * t263 + t268 * t297, -t244 * r_i_i_C(1) - (t266 * t261 - t264 * t300) * qJD(4) + t277 * t248 + t281 * t247, t247, -t262 * t297 - t263 * t295, 0; 0, t268 * qJD(2) - t306 * t265 - t290, -t303 * t295 + (t271 * qJD(3) + t289) * t262, t261 * t295 + t262 * t293, -t296, 0;];
JaD_transl  = t1;
