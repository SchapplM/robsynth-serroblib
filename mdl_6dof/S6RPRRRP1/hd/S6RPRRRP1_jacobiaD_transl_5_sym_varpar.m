% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:08:02
% EndTime: 2019-02-26 21:08:02
% DurationCPUTime: 0.24s
% Computational Cost: add. (382->56), mult. (398->81), div. (0->0), fcn. (299->10), ass. (0->53)
t265 = qJ(3) + qJ(4);
t261 = sin(t265);
t268 = cos(qJ(5));
t315 = r_i_i_C(1) * t268 + pkin(4);
t280 = t315 * t261;
t262 = cos(t265);
t298 = qJD(5) * t268;
t263 = qJD(3) + qJD(4);
t266 = sin(qJ(5));
t304 = t263 * t266;
t318 = t261 * t298 + t262 * t304;
t310 = pkin(9) + r_i_i_C(3);
t294 = t310 * t262;
t317 = (-pkin(4) * t261 + t294) * t263;
t267 = sin(qJ(3));
t306 = pkin(3) * qJD(3);
t296 = t267 * t306;
t316 = -t296 + t317;
t299 = qJD(5) * t266;
t287 = t261 * t299;
t313 = r_i_i_C(1) * t287 + t318 * r_i_i_C(2);
t300 = qJD(1) * t262;
t281 = -qJD(5) + t300;
t312 = t268 * t281;
t282 = qJD(5) * t262 - qJD(1);
t293 = t261 * t304;
t311 = t282 * t268 - t293;
t309 = pkin(3) * t267;
t303 = t263 * t268;
t264 = qJ(1) + pkin(10);
t259 = sin(t264);
t302 = qJD(1) * t259;
t260 = cos(t264);
t301 = qJD(1) * t260;
t297 = r_i_i_C(2) * t261 * t266;
t295 = t310 * t261;
t292 = t261 * t303;
t279 = t313 * t260 + t302 * t280;
t278 = t281 * t266;
t277 = t310 * t260 * t300 + t313 * t259 + t297 * t301;
t276 = -t294 - t297;
t269 = cos(qJ(3));
t275 = -t269 * pkin(3) - pkin(4) * t262 - pkin(2) - t295;
t274 = t282 * t266 + t292;
t273 = (-t262 * t315 - t295) * t263;
t272 = -t269 * t306 + t273;
t271 = (-t262 * t299 - t292) * r_i_i_C(1) + (-t262 * t298 + t293) * r_i_i_C(2) + t317;
t270 = -pkin(8) - pkin(7);
t243 = t274 * t259 - t260 * t312;
t242 = t311 * t259 + t260 * t278;
t241 = t259 * t312 + t274 * t260;
t240 = t259 * t278 - t311 * t260;
t1 = [t243 * r_i_i_C(1) + t242 * r_i_i_C(2) - t316 * t259 + (-cos(qJ(1)) * pkin(1) + t259 * t270 + t275 * t260) * qJD(1), 0 (t276 + t309) * t302 + t272 * t260 + t279, t260 * t273 + t276 * t302 + t279, t240 * r_i_i_C(1) + t241 * r_i_i_C(2), 0; -t241 * r_i_i_C(1) + t240 * r_i_i_C(2) + t316 * t260 + (-sin(qJ(1)) * pkin(1) - t260 * t270 + t275 * t259) * qJD(1), 0 (-t280 - t309) * t301 + t272 * t259 + t277, t259 * t273 - t280 * t301 + t277, -t242 * r_i_i_C(1) + t243 * r_i_i_C(2), 0; 0, 0, t271 - t296, t271 (-t262 * t303 + t287) * r_i_i_C(2) - t318 * r_i_i_C(1), 0;];
JaD_transl  = t1;
