% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:13
	% EndTime: 2019-12-29 20:05:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:18
	% EndTime: 2019-12-29 20:05:18
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:13
	% EndTime: 2019-12-29 20:05:13
	% DurationCPUTime: 0.14s
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
	t1 = [t23 * t25 + (-t18 * t28 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t22 * t18 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:13
	% EndTime: 2019-12-29 20:05:13
	% DurationCPUTime: 0.20s
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
	t1 = [t63 * t41 + (-t58 * t41 + t48 * t43) * qJD(1), (t40 * t56 - t43 * t54) * pkin(2) + t46, t46, 0, 0; -t63 * t43 + (t48 * t41 + t58 * t43) * qJD(1), t31 + (-pkin(2) * t54 - t53) * t41 + (-t49 - t62) * t55, -t43 * t51 + t31 + (-t36 * t55 - t41 * t57) * r_i_i_C(1), 0, 0; 0, -t63, -t47, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:15
	% EndTime: 2019-12-29 20:05:15
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (217->43), mult. (285->60), div. (0->0), fcn. (208->8), ass. (0->38)
	t224 = qJ(2) + qJ(3);
	t222 = cos(t224);
	t259 = r_i_i_C(3) + qJ(4);
	t270 = t259 * t222;
	t221 = sin(t224);
	t219 = t221 * qJD(4);
	t223 = qJD(2) + qJD(3);
	t225 = sin(pkin(9));
	t226 = cos(pkin(9));
	t260 = r_i_i_C(2) * t225;
	t268 = r_i_i_C(1) * t226 + pkin(3);
	t237 = t268 - t260;
	t227 = sin(qJ(2));
	t258 = pkin(2) * qJD(2);
	t250 = t227 * t258;
	t269 = (-t237 * t221 + t270) * t223 + (t225 * r_i_i_C(1) + t226 * r_i_i_C(2) + pkin(6) + pkin(7)) * qJD(1) + t219 - t250;
	t251 = t223 * t260;
	t267 = (qJD(4) + t251) * t222;
	t262 = pkin(2) * t227;
	t228 = sin(qJ(1));
	t256 = t222 * t228;
	t230 = cos(qJ(1));
	t255 = t223 * t230;
	t254 = qJD(1) * t228;
	t253 = qJD(1) * t230;
	t248 = t221 * t254;
	t247 = t221 * t253;
	t245 = t259 * t221;
	t243 = t259 * t228;
	t241 = t267 * t230 + t268 * t248;
	t240 = t268 * t223;
	t239 = t268 * t230;
	t238 = t267 * t228 + t247 * t260 + t253 * t270;
	t234 = t219 + t223 * t270 + (-t240 + t251) * t221;
	t229 = cos(qJ(2));
	t233 = qJD(1) * (-t229 * pkin(2) - t237 * t222 - pkin(1) - t245);
	t232 = -t229 * t258 + (-t222 * t268 - t245) * t223;
	t1 = [-t269 * t228 + t230 * t233, (-t221 * t260 + t262 - t270) * t254 + t232 * t230 + t241, (-t254 * t260 - t259 * t255) * t221 + (-qJD(1) * t243 - t223 * t239) * t222 + t241, t222 * t255 - t248, 0; t228 * t233 + t269 * t230, (-t221 * t268 - t262) * t253 + t232 * t228 + t238, -t240 * t256 + (-qJD(1) * t239 - t223 * t243) * t221 + t238, t223 * t256 + t247, 0; 0, t234 - t250, t234, t223 * t221, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:15
	% EndTime: 2019-12-29 20:05:16
	% DurationCPUTime: 0.75s
	% Computational Cost: add. (402->68), mult. (437->101), div. (0->0), fcn. (338->10), ass. (0->56)
	t262 = pkin(9) + qJ(5);
	t258 = sin(t262);
	t259 = cos(t262);
	t264 = qJ(2) + qJ(3);
	t260 = sin(t264);
	t299 = qJD(5) * t260;
	t261 = cos(t264);
	t263 = qJD(2) + qJD(3);
	t306 = t261 * t263;
	t321 = t258 * t306 + t259 * t299;
	t256 = cos(pkin(9)) * pkin(4) + pkin(3);
	t320 = r_i_i_C(1) * t259 + t256;
	t255 = t260 * qJD(4);
	t267 = sin(qJ(2));
	t308 = pkin(2) * qJD(2);
	t297 = t267 * t308;
	t266 = -pkin(8) - qJ(4);
	t309 = r_i_i_C(3) - t266;
	t319 = (-t256 * t260 + t261 * t309) * t263 + (pkin(4) * sin(pkin(9)) + pkin(7) + pkin(6)) * qJD(1) + t255 - t297;
	t291 = t258 * t299;
	t318 = r_i_i_C(1) * t291 + r_i_i_C(2) * t321 + qJD(4) * t261;
	t270 = cos(qJ(1));
	t284 = qJD(5) * t261 - qJD(1);
	t317 = t270 * t284;
	t283 = qJD(1) * t261 - qJD(5);
	t268 = sin(qJ(1));
	t304 = t263 * t268;
	t295 = t260 * t304;
	t315 = t270 * t283 - t295;
	t313 = pkin(2) * t267;
	t311 = r_i_i_C(2) * t258;
	t310 = r_i_i_C(3) * t261;
	t305 = t261 * t266;
	t303 = t263 * t270;
	t302 = qJD(1) * t268;
	t301 = qJD(1) * t270;
	t298 = t260 * t311;
	t294 = t260 * t303;
	t293 = t260 * t302;
	t292 = t260 * t301;
	t282 = t320 * t263;
	t281 = t284 * t268;
	t279 = t266 * t295 + t268 * t318 + t292 * t311 + t301 * t310;
	t278 = -t260 * t320 - t305;
	t277 = t266 * t294 + t270 * t318 + t293 * t320 + t302 * t305;
	t269 = cos(qJ(2));
	t276 = qJD(1) * (-pkin(2) * t269 - t256 * t261 - t260 * t309 - pkin(1));
	t275 = (-r_i_i_C(3) * t260 - t261 * t320) * t263;
	t274 = t268 * t283 + t294;
	t273 = -t269 * t308 + t275;
	t272 = t263 * t298 + r_i_i_C(3) * t306 + t255 - t260 * t282 + (-t263 * t266 + (-r_i_i_C(1) * t258 - r_i_i_C(2) * t259) * qJD(5)) * t261;
	t237 = t258 * t281 - t259 * t315;
	t236 = t258 * t315 + t259 * t281;
	t235 = t258 * t317 + t259 * t274;
	t234 = t258 * t274 - t259 * t317;
	t1 = [t237 * r_i_i_C(1) + t236 * r_i_i_C(2) - t268 * t319 + t270 * t276, (-t298 - t310 + t313) * t302 + t273 * t270 + t277, (-r_i_i_C(3) * t303 - t302 * t311) * t260 + (-r_i_i_C(3) * t302 - t270 * t282) * t261 + t277, t261 * t303 - t293, r_i_i_C(1) * t234 + r_i_i_C(2) * t235; -t235 * r_i_i_C(1) + t234 * r_i_i_C(2) + t268 * t276 + t270 * t319, t273 * t268 + (t278 - t313) * t301 + t279, t268 * t275 + t278 * t301 + t279, t261 * t304 + t292, -r_i_i_C(1) * t236 + r_i_i_C(2) * t237; 0, t272 - t297, t272, t263 * t260, (-t259 * t306 + t291) * r_i_i_C(2) - t321 * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end