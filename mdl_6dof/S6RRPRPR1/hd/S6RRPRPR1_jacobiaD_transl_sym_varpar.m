% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:04
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(7) + r_i_i_C(3);
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
	t1 = [t23 * t25 + (-t28 * t18 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0, 0; -t20 * t21 + (t22 * t18 + t28 * t20) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0, 0; 0, -t21, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (50->16), mult. (90->25), div. (0->0), fcn. (59->6), ass. (0->14)
	t27 = qJ(2) + pkin(10);
	t25 = sin(t27);
	t26 = cos(t27);
	t42 = -r_i_i_C(1) * t26 + r_i_i_C(2) * t25 - cos(qJ(2)) * pkin(2);
	t40 = r_i_i_C(3) + qJ(3) + pkin(7);
	t32 = cos(qJ(1));
	t39 = qJD(1) * t32;
	t37 = -pkin(1) + t42;
	t36 = sin(qJ(2)) * pkin(2) + r_i_i_C(1) * t25 + r_i_i_C(2) * t26;
	t30 = sin(qJ(1));
	t35 = t36 * t30;
	t34 = qJD(2) * t42;
	t33 = t36 * qJD(2);
	t1 = [t32 * qJD(3) + qJD(2) * t35 + (-t40 * t30 + t37 * t32) * qJD(1), qJD(1) * t35 + t32 * t34, t39, 0, 0, 0; t30 * qJD(3) - t32 * t33 + (t37 * t30 + t40 * t32) * qJD(1), t30 * t34 - t36 * t39, qJD(1) * t30, 0, 0, 0; 0, -t33, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (131->30), mult. (132->40), div. (0->0), fcn. (86->8), ass. (0->27)
	t49 = qJ(2) + pkin(10);
	t41 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t49);
	t48 = qJD(2) + qJD(4);
	t46 = qJ(4) + t49;
	t43 = cos(t46);
	t67 = r_i_i_C(2) * t43;
	t42 = sin(t46);
	t69 = r_i_i_C(1) * t42;
	t56 = t67 + t69;
	t54 = t56 * t48;
	t70 = t41 * qJD(2) - t54;
	t68 = r_i_i_C(2) * t42;
	t66 = r_i_i_C(3) + pkin(8) + qJ(3) + pkin(7);
	t65 = t43 * t48;
	t51 = sin(qJ(1));
	t64 = qJD(1) * t51;
	t53 = cos(qJ(1));
	t63 = qJD(1) * t53;
	t62 = r_i_i_C(1) * t65;
	t61 = t48 * t68;
	t59 = qJD(1) * t67;
	t60 = t51 * t59 + t53 * t61 + t64 * t69;
	t57 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t49);
	t58 = t57 * qJD(2) - t62;
	t55 = -r_i_i_C(1) * t43 - pkin(1) + t57 + t68;
	t36 = t51 * t61;
	t1 = [t53 * qJD(3) - t70 * t51 + (-t66 * t51 + t55 * t53) * qJD(1), -t41 * t64 + t58 * t53 + t60, t63, -t53 * t62 + t60, 0, 0; t51 * qJD(3) + t70 * t53 + (t55 * t51 + t66 * t53) * qJD(1), t36 + t58 * t51 + (t41 - t56) * t63, t64, -t53 * t59 + t36 + (-t42 * t63 - t51 * t65) * r_i_i_C(1), 0, 0; 0, t70, 0, -t54, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:27
	% EndTime: 2019-10-10 10:04:27
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (322->48), mult. (303->61), div. (0->0), fcn. (221->10), ass. (0->37)
	t234 = qJ(2) + pkin(10);
	t231 = qJ(4) + t234;
	t228 = cos(t231);
	t267 = r_i_i_C(3) + qJ(5);
	t277 = t267 * t228;
	t222 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t234);
	t214 = t222 * qJD(2);
	t227 = sin(t231);
	t226 = t227 * qJD(5);
	t233 = qJD(2) + qJD(4);
	t235 = sin(pkin(11));
	t236 = cos(pkin(11));
	t268 = r_i_i_C(2) * t235;
	t275 = r_i_i_C(1) * t236 + pkin(4);
	t245 = t275 - t268;
	t276 = (-t245 * t227 + t277) * t233 + (t235 * r_i_i_C(1) + t236 * r_i_i_C(2) + pkin(7) + pkin(8) + qJ(3)) * qJD(1) + t214 + t226;
	t260 = t233 * t268;
	t274 = (qJD(5) + t260) * t228;
	t238 = sin(qJ(1));
	t265 = t228 * t238;
	t240 = cos(qJ(1));
	t264 = t233 * t240;
	t263 = qJD(1) * t238;
	t262 = qJD(1) * t240;
	t258 = t227 * t263;
	t257 = t227 * t262;
	t255 = t267 * t227;
	t253 = t267 * t238;
	t250 = t274 * t240 + t275 * t258;
	t249 = t275 * t233;
	t248 = t275 * t240;
	t247 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t234);
	t246 = t274 * t238 + t257 * t268 + t262 * t277;
	t243 = t226 + t233 * t277 + (-t249 + t260) * t227;
	t242 = t247 * qJD(2) + (-t228 * t275 - t255) * t233;
	t241 = qJD(3) + (-t245 * t228 - pkin(1) + t247 - t255) * qJD(1);
	t1 = [-t276 * t238 + t241 * t240, (-t227 * t268 - t222 - t277) * t263 + t242 * t240 + t250, t262, (-t263 * t268 - t267 * t264) * t227 + (-qJD(1) * t253 - t233 * t248) * t228 + t250, t228 * t264 - t258, 0; t241 * t238 + t276 * t240, (-t227 * t275 + t222) * t262 + t242 * t238 + t246, t263, -t249 * t265 + (-qJD(1) * t248 - t233 * t253) * t227 + t246, t233 * t265 + t257, 0; 0, t214 + t243, 0, t243, t233 * t227, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:27
	% EndTime: 2019-10-10 10:04:27
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (539->73), mult. (455->102), div. (0->0), fcn. (351->12), ass. (0->55)
	t276 = pkin(11) + qJ(6);
	t270 = sin(t276);
	t272 = cos(t276);
	t278 = qJ(2) + pkin(10);
	t274 = qJ(4) + t278;
	t267 = sin(t274);
	t312 = qJD(6) * t267;
	t268 = cos(t274);
	t277 = qJD(2) + qJD(4);
	t319 = t268 * t277;
	t332 = t270 * t319 + t272 * t312;
	t269 = cos(pkin(11)) * pkin(5) + pkin(4);
	t331 = r_i_i_C(1) * t272 + t269;
	t262 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t278);
	t257 = t262 * qJD(2);
	t266 = t267 * qJD(5);
	t280 = -pkin(9) - qJ(5);
	t321 = r_i_i_C(3) - t280;
	t330 = (-t267 * t269 + t321 * t268) * t277 + (pkin(5) * sin(pkin(11)) + pkin(7) + pkin(8) + qJ(3)) * qJD(1) + t257 + t266;
	t305 = t270 * t312;
	t329 = r_i_i_C(1) * t305 + t332 * r_i_i_C(2) + qJD(5) * t268;
	t284 = cos(qJ(1));
	t297 = qJD(6) * t268 - qJD(1);
	t328 = t284 * t297;
	t296 = qJD(1) * t268 - qJD(6);
	t282 = sin(qJ(1));
	t317 = t277 * t282;
	t310 = t267 * t317;
	t326 = t296 * t284 - t310;
	t323 = r_i_i_C(2) * t270;
	t322 = r_i_i_C(3) * t268;
	t318 = t268 * t280;
	t316 = t277 * t284;
	t315 = qJD(1) * t282;
	t314 = qJD(1) * t284;
	t311 = t267 * t323;
	t309 = t267 * t316;
	t307 = t267 * t315;
	t306 = t267 * t314;
	t295 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t278);
	t294 = t331 * t277;
	t293 = t297 * t282;
	t292 = t280 * t310 + t329 * t282 + t306 * t323 + t314 * t322;
	t291 = -t267 * t331 - t318;
	t290 = t280 * t309 + t329 * t284 + t331 * t307 + t315 * t318;
	t289 = (-r_i_i_C(3) * t267 - t268 * t331) * t277;
	t288 = t296 * t282 + t309;
	t287 = qJD(3) + (-t321 * t267 - t268 * t269 - pkin(1) + t295) * qJD(1);
	t286 = t295 * qJD(2) + t289;
	t285 = t277 * t311 + r_i_i_C(3) * t319 + t266 - t267 * t294 + (-t277 * t280 + (-r_i_i_C(1) * t270 - r_i_i_C(2) * t272) * qJD(6)) * t268;
	t244 = t270 * t293 - t326 * t272;
	t243 = t326 * t270 + t272 * t293;
	t242 = t270 * t328 + t288 * t272;
	t241 = t288 * t270 - t272 * t328;
	t1 = [t244 * r_i_i_C(1) + t243 * r_i_i_C(2) - t330 * t282 + t287 * t284, (-t262 - t311 - t322) * t315 + t286 * t284 + t290, t314, (-r_i_i_C(3) * t316 - t315 * t323) * t267 + (-r_i_i_C(3) * t315 - t284 * t294) * t268 + t290, t268 * t316 - t307, t241 * r_i_i_C(1) + t242 * r_i_i_C(2); -t242 * r_i_i_C(1) + t241 * r_i_i_C(2) + t287 * t282 + t330 * t284, t286 * t282 + (t262 + t291) * t314 + t292, t315, t282 * t289 + t291 * t314 + t292, t268 * t317 + t306, -t243 * r_i_i_C(1) + t244 * r_i_i_C(2); 0, t257 + t285, 0, t285, t277 * t267, (-t272 * t319 + t305) * r_i_i_C(2) - t332 * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end