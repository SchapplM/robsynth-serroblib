% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRP6
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRP6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRP6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:55:14
	% EndTime: 2019-12-31 21:55:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:55:14
	% EndTime: 2019-12-31 21:55:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:55:14
	% EndTime: 2019-12-31 21:55:14
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(6) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t18 * t28 + t20 * t22) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t18 * t22 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:55:14
	% EndTime: 2019-12-31 21:55:14
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (81->26), mult. (114->37), div. (0->0), fcn. (73->6), ass. (0->27)
	t38 = qJD(2) + qJD(3);
	t39 = qJ(2) + qJ(3);
	t37 = cos(t39);
	t59 = r_i_i_C(2) * t37;
	t36 = sin(t39);
	t61 = r_i_i_C(1) * t36;
	t49 = t59 + t61;
	t47 = t49 * t38;
	t40 = sin(qJ(2));
	t62 = pkin(2) * t40;
	t63 = qJD(2) * t62 + t47;
	t60 = r_i_i_C(2) * t36;
	t58 = r_i_i_C(3) + pkin(7) + pkin(6);
	t57 = t37 * t38;
	t41 = sin(qJ(1));
	t56 = qJD(1) * t41;
	t43 = cos(qJ(1));
	t55 = qJD(1) * t43;
	t42 = cos(qJ(2));
	t54 = qJD(2) * t42;
	t53 = r_i_i_C(1) * t57;
	t52 = t38 * t60;
	t51 = qJD(1) * t59;
	t48 = -t42 * pkin(2) - r_i_i_C(1) * t37 - pkin(1) + t60;
	t46 = t41 * t51 + t56 * t61 + (t52 - t53) * t43;
	t31 = t41 * t52;
	t1 = [t63 * t41 + (-t41 * t58 + t43 * t48) * qJD(1), (t40 * t56 - t43 * t54) * pkin(2) + t46, t46, 0, 0; -t63 * t43 + (t41 * t48 + t43 * t58) * qJD(1), t31 + (-pkin(2) * t54 - t53) * t41 + (-t49 - t62) * t55, -t43 * t51 + t31 + (-t36 * t55 - t41 * t57) * r_i_i_C(1), 0, 0; 0, -t63, -t47, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:55:15
	% EndTime: 2019-12-31 21:55:16
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (270->58), mult. (394->95), div. (0->0), fcn. (297->8), ass. (0->54)
	t258 = qJ(2) + qJ(3);
	t255 = sin(t258);
	t262 = cos(qJ(4));
	t310 = r_i_i_C(1) * t262 + pkin(3);
	t313 = t255 * t310;
	t259 = sin(qJ(4));
	t293 = qJD(4) * t262;
	t256 = cos(t258);
	t257 = qJD(2) + qJD(3);
	t301 = t256 * t257;
	t312 = t255 * t293 + t259 * t301;
	t307 = pkin(8) + r_i_i_C(3);
	t288 = t307 * t256;
	t260 = sin(qJ(2));
	t302 = pkin(2) * qJD(2);
	t290 = t260 * t302;
	t305 = pkin(3) * t255;
	t311 = -t290 + (t288 - t305) * t257;
	t294 = qJD(4) * t259;
	t281 = t255 * t294;
	t308 = r_i_i_C(1) * t281 + r_i_i_C(2) * t312;
	t306 = pkin(2) * t260;
	t303 = r_i_i_C(2) * t259;
	t261 = sin(qJ(1));
	t300 = t257 * t261;
	t299 = t257 * t262;
	t264 = cos(qJ(1));
	t298 = t257 * t264;
	t297 = t262 * t264;
	t296 = qJD(1) * t261;
	t295 = qJD(1) * t264;
	t292 = t255 * t303;
	t291 = qJD(1) * t303;
	t289 = t307 * t255;
	t287 = t307 * t261;
	t286 = t255 * t299;
	t276 = qJD(4) * t256 - qJD(1);
	t275 = qJD(1) * t256 - qJD(4);
	t274 = t310 * t256;
	t273 = t310 * t264;
	t272 = t264 * t308 + t296 * t313;
	t271 = t276 * t259;
	t270 = t264 * t255 * t291 + t261 * t308 + t288 * t295;
	t263 = cos(qJ(2));
	t269 = -pkin(2) * t263 - pkin(3) * t256 - pkin(1) - t289;
	t268 = t255 * t298 + t261 * t275;
	t267 = -t263 * t302 + (-t274 - t289) * t257;
	t266 = -t256 * r_i_i_C(2) * t293 + (-t256 * t294 - t286) * r_i_i_C(1) + t307 * t301 + (-t305 + t292) * t257;
	t265 = -pkin(7) - pkin(6);
	t239 = -t275 * t297 + (t271 + t286) * t261;
	t238 = t276 * t262 * t261 + (-t255 * t300 + t264 * t275) * t259;
	t237 = t262 * t268 + t264 * t271;
	t236 = t259 * t268 - t276 * t297;
	t1 = [t239 * r_i_i_C(1) + t238 * r_i_i_C(2) - t311 * t261 + (t261 * t265 + t264 * t269) * qJD(1), (-t288 - t292 + t306) * t296 + t267 * t264 + t272, (-t261 * t291 - t298 * t307) * t255 + (-qJD(1) * t287 - t257 * t273) * t256 + t272, r_i_i_C(1) * t236 + r_i_i_C(2) * t237, 0; -t237 * r_i_i_C(1) + t236 * r_i_i_C(2) + t311 * t264 + (t261 * t269 - t264 * t265) * qJD(1), (-t306 - t313) * t295 + t267 * t261 + t270, -t274 * t300 + (-qJD(1) * t273 - t257 * t287) * t255 + t270, -r_i_i_C(1) * t238 + r_i_i_C(2) * t239, 0; 0, t266 - t290, t266, (-t256 * t299 + t281) * r_i_i_C(2) - t312 * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:55:15
	% EndTime: 2019-12-31 21:55:16
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (367->70), mult. (515->104), div. (0->0), fcn. (394->8), ass. (0->56)
	t322 = pkin(4) + r_i_i_C(1);
	t262 = qJD(2) + qJD(3);
	t265 = sin(qJ(4));
	t315 = r_i_i_C(2) * t265;
	t329 = t262 * t315 + qJD(5);
	t263 = qJ(2) + qJ(3);
	t260 = sin(t263);
	t261 = cos(t263);
	t268 = cos(qJ(4));
	t303 = qJD(4) * t268;
	t304 = qJD(4) * t265;
	t328 = t329 * t261 + (r_i_i_C(2) * t303 + t304 * t322) * t260;
	t267 = sin(qJ(1));
	t270 = cos(qJ(1));
	t287 = qJD(1) * t261 - qJD(4);
	t281 = t287 * t270;
	t288 = qJD(4) * t261 - qJD(1);
	t283 = t288 * t268;
	t308 = t262 * t267;
	t298 = t260 * t308;
	t236 = t267 * t283 + (t281 - t298) * t265;
	t258 = t268 * pkin(4) + pkin(3);
	t327 = r_i_i_C(1) * t268 + t258;
	t266 = sin(qJ(2));
	t312 = pkin(2) * qJD(2);
	t299 = t266 * t312;
	t264 = -qJ(5) - pkin(8);
	t313 = r_i_i_C(3) - t264;
	t323 = (-pkin(4) * t304 + t313 * t262) * t261 + (pkin(4) * t265 + pkin(6) + pkin(7)) * qJD(1) - (t258 * t262 - qJD(5)) * t260 - t299;
	t317 = pkin(2) * t266;
	t314 = r_i_i_C(3) * t261;
	t311 = t261 * t262;
	t310 = t261 * t264;
	t309 = t262 * t260;
	t307 = t262 * t270;
	t306 = qJD(1) * t267;
	t305 = qJD(1) * t270;
	t297 = t260 * t307;
	t296 = t260 * t306;
	t295 = t260 * t305;
	t284 = t327 * t262;
	t282 = t288 * t265;
	t280 = -r_i_i_C(2) * t268 - t265 * t322;
	t279 = -t260 * t327 - t310;
	t278 = t264 * t298 + t328 * t267 + t295 * t315 + t305 * t314;
	t277 = t264 * t297 + t328 * t270 + t327 * t296 + t306 * t310;
	t276 = (-r_i_i_C(3) * t260 - t261 * t327) * t262;
	t275 = t287 * t267 + t297;
	t269 = cos(qJ(2));
	t274 = pkin(4) * t303 + (-t269 * pkin(2) - t258 * t261 - t313 * t260 - pkin(1)) * qJD(1);
	t273 = -t269 * t312 + t276;
	t234 = t275 * t265 - t270 * t283;
	t272 = r_i_i_C(3) * t311 + (t280 * qJD(4) - t262 * t264) * t261 + (-t284 + t329) * t260;
	t237 = -t268 * t281 + (t268 * t309 + t282) * t267;
	t235 = t275 * t268 + t270 * t282;
	t1 = [t237 * r_i_i_C(1) + t236 * r_i_i_C(2) - t323 * t267 + t274 * t270, (-t260 * t315 - t314 + t317) * t306 + t273 * t270 + t277, (-r_i_i_C(3) * t307 - t306 * t315) * t260 + (-r_i_i_C(3) * t306 - t270 * t284) * t261 + t277, t235 * r_i_i_C(2) + t322 * t234, t261 * t307 - t296; -t235 * r_i_i_C(1) + t234 * r_i_i_C(2) + t274 * t267 + t323 * t270, t273 * t267 + (t279 - t317) * t305 + t278, t267 * t276 + t279 * t305 + t278, t237 * r_i_i_C(2) - t322 * t236, t261 * t308 + t295; 0, t272 - t299, t272, t280 * t311 + (-t268 * t322 + t315) * t260 * qJD(4), t309;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end