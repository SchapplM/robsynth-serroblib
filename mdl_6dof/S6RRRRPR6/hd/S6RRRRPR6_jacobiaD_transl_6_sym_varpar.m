% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:33:31
% EndTime: 2019-02-26 22:33:31
% DurationCPUTime: 0.32s
% Computational Cost: add. (845->75), mult. (649->100), div. (0->0), fcn. (512->12), ass. (0->63)
t260 = qJ(3) + qJ(4);
t254 = pkin(11) + t260;
t244 = -pkin(4) * sin(t260) - pkin(5) * sin(t254);
t259 = qJD(3) + qJD(4);
t238 = t244 * t259;
t261 = sin(qJ(3));
t298 = pkin(3) * qJD(3);
t236 = -t261 * t298 + t238;
t264 = cos(qJ(3));
t307 = pkin(5) * cos(t254) + pkin(4) * cos(t260);
t242 = t264 * pkin(3) + t307;
t240 = pkin(2) + t242;
t241 = t261 * pkin(3) - t244;
t262 = sin(qJ(2));
t265 = cos(qJ(2));
t287 = t262 * qJD(5);
t299 = r_i_i_C(3) + pkin(10) + qJ(5) + pkin(9) + pkin(8);
t308 = t299 * t265;
t312 = (-t240 * t262 + t308) * qJD(2) + (pkin(7) + t241) * qJD(1) + t265 * t236 + t287;
t252 = qJ(6) + t254;
t247 = sin(t252);
t302 = r_i_i_C(2) * t247;
t248 = cos(t252);
t303 = r_i_i_C(1) * t248;
t273 = t240 - t302 + t303;
t268 = -t273 * t262 + t308;
t266 = cos(qJ(1));
t253 = qJD(6) + t259;
t281 = t253 * t265 - qJD(1);
t310 = t266 * t281;
t239 = t307 * t259;
t301 = r_i_i_C(2) * t248;
t278 = r_i_i_C(1) * t247 + t301;
t306 = t253 * t278 - t236;
t263 = sin(qJ(1));
t292 = qJD(1) * t265;
t280 = -t253 + t292;
t290 = qJD(2) * t262;
t305 = -t263 * t290 + t280 * t266;
t296 = t253 * t262;
t288 = qJD(2) * t266;
t270 = t262 * t288 + t280 * t263;
t232 = t270 * t247 - t248 * t310;
t233 = t247 * t310 + t270 * t248;
t295 = t232 * r_i_i_C(1) + t233 * r_i_i_C(2);
t275 = t281 * t263;
t234 = t305 * t247 + t248 * t275;
t235 = t247 * t275 - t305 * t248;
t294 = -t234 * r_i_i_C(1) + t235 * r_i_i_C(2);
t293 = qJD(1) * t263;
t291 = qJD(1) * t266;
t289 = qJD(2) * t265;
t286 = t253 * t303;
t284 = t299 * t262;
t277 = t241 * t292 + t236;
t276 = t244 * t292 - t238;
t237 = t264 * t298 + t239;
t272 = qJD(1) * t242 - t237 * t265 + t241 * t290;
t271 = qJD(1) * t307 - t239 * t265 - t244 * t290;
t269 = t237 + (-t240 * t265 - pkin(1) - t284) * qJD(1);
t267 = qJD(5) * t265 + t306 * t262 + (-t273 * t265 - t284) * qJD(2);
t243 = t296 * t302;
t1 = [t235 * r_i_i_C(1) + t234 * r_i_i_C(2) - t312 * t263 + t269 * t266, t267 * t266 - t268 * t293, t263 * t277 + t266 * t272 + t295, -t263 * t276 + t266 * t271 + t295, -t262 * t293 + t265 * t288, t295; -t233 * r_i_i_C(1) + t232 * r_i_i_C(2) + t269 * t263 + t312 * t266, t263 * t267 + t268 * t291, t263 * t272 - t266 * t277 + t294, t263 * t271 + t266 * t276 + t294, t262 * t291 + t263 * t289, t294; 0, t268 * qJD(2) - t306 * t265 + t287, t243 + (-t237 - t286) * t262 + (-t241 - t278) * t289, t243 + (-t239 - t286) * t262 + (t244 - t278) * t289, t290, -t289 * t301 + t243 + (-t247 * t289 - t248 * t296) * r_i_i_C(1);];
JaD_transl  = t1;
