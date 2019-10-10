% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:33
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
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
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) + r_i_i_C(2);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (23->16), mult. (72->29), div. (0->0), fcn. (46->4), ass. (0->13)
	t17 = sin(qJ(3));
	t19 = cos(qJ(3));
	t29 = (r_i_i_C(1) * t19 - r_i_i_C(2) * t17) * qJD(3);
	t18 = sin(qJ(1));
	t28 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t27 = qJD(1) * t20;
	t26 = qJD(3) * t18;
	t25 = qJD(3) * t20;
	t24 = -pkin(1) - pkin(7) - r_i_i_C(3);
	t22 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19 + qJ(2);
	t21 = qJD(2) + t29;
	t1 = [t21 * t20 + (-t18 * t22 + t20 * t24) * qJD(1), t27, (-t17 * t27 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 + t19 * t27) * r_i_i_C(1), 0, 0, 0; t21 * t18 + (t18 * t24 + t20 * t22) * qJD(1), t28, (-t17 * t28 + t19 * t25) * r_i_i_C(2) + (t17 * t25 + t19 * t28) * r_i_i_C(1), 0, 0, 0; 0, 0, -t29, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (85->27), mult. (126->41), div. (0->0), fcn. (81->6), ass. (0->27)
	t43 = cos(qJ(3));
	t62 = pkin(3) * t43;
	t40 = qJ(3) + qJ(4);
	t38 = cos(t40);
	t61 = r_i_i_C(1) * t38;
	t37 = sin(t40);
	t60 = r_i_i_C(2) * t37;
	t39 = qJD(3) + qJD(4);
	t59 = t37 * t39;
	t58 = t38 * t39;
	t42 = sin(qJ(1));
	t57 = qJD(1) * t42;
	t44 = cos(qJ(1));
	t56 = qJD(1) * t44;
	t41 = sin(qJ(3));
	t55 = qJD(3) * t41;
	t54 = -pkin(1) - r_i_i_C(3) - pkin(8) - pkin(7);
	t53 = r_i_i_C(1) * t59;
	t52 = qJD(1) * t61;
	t51 = qJD(3) * t62;
	t50 = -r_i_i_C(1) * t58 + r_i_i_C(2) * t59;
	t49 = -r_i_i_C(1) * t37 - r_i_i_C(2) * t38;
	t48 = t42 * t52 - t57 * t60 + (t58 * r_i_i_C(2) + t53) * t44;
	t47 = pkin(3) * t41 + qJ(2) - t49;
	t46 = t51 + qJD(2) + (-t60 + t61) * t39;
	t35 = t44 * t52;
	t1 = [t46 * t44 + (-t47 * t42 + t54 * t44) * qJD(1), t56, t35 + (-t60 + t62) * t56 + (-pkin(3) * t55 + t49 * t39) * t42, -t42 * t53 + t35 + (-t37 * t56 - t42 * t58) * r_i_i_C(2), 0, 0; t46 * t42 + (t54 * t42 + t47 * t44) * qJD(1), t57, (t43 * t57 + t44 * t55) * pkin(3) + t48, t48, 0, 0; 0, 0, t50 - t51, t50, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (154->35), mult. (158->46), div. (0->0), fcn. (103->8), ass. (0->32)
	t56 = qJ(3) + qJ(4);
	t50 = pkin(10) + t56;
	t48 = sin(t50);
	t49 = cos(t50);
	t82 = r_i_i_C(1) * t48 + r_i_i_C(2) * t49;
	t53 = cos(t56);
	t74 = r_i_i_C(2) * t48;
	t80 = pkin(4) * t53 - t74;
	t52 = sin(t56);
	t81 = -pkin(4) * t52 - t82;
	t55 = qJD(3) + qJD(4);
	t77 = pkin(4) * t55;
	t75 = r_i_i_C(1) * t49;
	t60 = cos(qJ(1));
	t72 = t55 * t60;
	t71 = pkin(3) * qJD(3);
	t58 = sin(qJ(1));
	t70 = qJD(1) * t58;
	t51 = qJD(1) * t60;
	t69 = -pkin(1) - r_i_i_C(3) - qJ(5) - pkin(8) - pkin(7);
	t67 = qJD(1) * t75;
	t68 = t58 * t67 + t82 * t72;
	t59 = cos(qJ(3));
	t66 = t59 * t71;
	t64 = qJD(1) * (pkin(3) * t59 + t80);
	t57 = sin(qJ(3));
	t63 = pkin(3) * t57 + qJ(2) - t81;
	t62 = (-t75 - t80) * t55;
	t61 = qJD(2) + t53 * t77 + t66 + (-t74 + t75) * t55;
	t44 = t60 * t67;
	t39 = -t52 * t77 - t57 * t71;
	t1 = [-t58 * qJD(5) + t61 * t60 + (-t63 * t58 + t69 * t60) * qJD(1), t51, t44 + t60 * t64 + (-t55 * t82 + t39) * t58, t81 * t58 * t55 + t80 * t51 + t44, -t70, 0; qJD(5) * t60 + t61 * t58 + (t69 * t58 + t63 * t60) * qJD(1), t70, -t60 * t39 + t58 * t64 + t68, -t70 * t74 + (t52 * t72 + t53 * t70) * pkin(4) + t68, t51, 0; 0, 0, t62 - t66, t62, 0, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:57
	% EndTime: 2019-10-10 01:33:57
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (417->71), mult. (438->105), div. (0->0), fcn. (327->10), ass. (0->58)
	t268 = qJ(3) + qJ(4);
	t262 = pkin(10) + t268;
	t261 = cos(t262);
	t272 = cos(qJ(6));
	t324 = r_i_i_C(1) * t272 + pkin(5);
	t325 = t261 * t324;
	t316 = pkin(9) + r_i_i_C(3);
	t323 = qJD(1) * t325;
	t270 = sin(qJ(3));
	t298 = t316 * t261;
	t260 = sin(t262);
	t312 = pkin(5) * t260;
	t264 = sin(t268);
	t315 = pkin(4) * t264;
	t321 = -t270 * pkin(3) - qJ(2) + t298 - t312 - t315;
	t265 = cos(t268);
	t269 = sin(qJ(6));
	t318 = -r_i_i_C(2) * t261 * t269 + pkin(4) * t265;
	t274 = cos(qJ(1));
	t285 = qJD(1) * t260 + qJD(6);
	t320 = t285 * t274;
	t319 = t316 * t260;
	t271 = sin(qJ(1));
	t267 = qJD(3) + qJD(4);
	t306 = t267 * t274;
	t317 = -t261 * t306 + t285 * t271;
	t313 = pkin(4) * t267;
	t310 = -pkin(1) - qJ(5) - pkin(8) - pkin(7);
	t309 = pkin(3) * qJD(3);
	t308 = t260 * t269;
	t307 = t267 * t272;
	t304 = qJD(1) * t271;
	t263 = qJD(1) * t274;
	t303 = qJD(6) * t269;
	t302 = qJD(6) * t272;
	t273 = cos(qJ(3));
	t299 = t273 * t309;
	t297 = t267 * t308;
	t295 = t261 * t267 * t271;
	t290 = t261 * t303;
	t289 = t261 * t302;
	t287 = r_i_i_C(2) * t289;
	t286 = -qJD(6) * t260 - qJD(1);
	t283 = t286 * t274;
	t282 = qJD(1) * (t273 * pkin(3) + t318);
	t281 = -r_i_i_C(2) * t308 - t298;
	t280 = t271 * r_i_i_C(2) * t297 + t323 * t274 + t316 * (t260 * t263 + t295);
	t279 = qJD(1) * t318;
	t278 = t260 * t307 + t290;
	t277 = (r_i_i_C(1) * t290 + t287) * t274 + t323 * t271 + (t316 * t304 + t324 * t306) * t260;
	t276 = qJD(2) + t265 * t313 + t299 + (pkin(5) * t261 + t319) * t267;
	t275 = (-t318 - t319 - t325) * t267 + (r_i_i_C(1) * t303 + r_i_i_C(2) * t302) * t260;
	t240 = -t264 * t313 - t270 * t309;
	t237 = t272 * t320 + (t261 * t307 + t286 * t269) * t271;
	t236 = t286 * t272 * t271 + (-t295 - t320) * t269;
	t235 = t269 * t283 - t317 * t272;
	t234 = t317 * t269 + t272 * t283;
	t1 = [t235 * r_i_i_C(1) + t234 * r_i_i_C(2) - t271 * qJD(5) + t276 * t274 + (t321 * t271 + t310 * t274) * qJD(1), t263, t274 * t282 + (-t278 * r_i_i_C(1) - t267 * t312 + t240 - t287) * t271 + t280, t274 * t279 + ((-r_i_i_C(1) * t269 - r_i_i_C(2) * t272) * t261 * qJD(6) + (-t260 * t324 - t315) * t267) * t271 + t280, -t304, t236 * r_i_i_C(1) - t237 * r_i_i_C(2); t237 * r_i_i_C(1) + t236 * r_i_i_C(2) + t274 * qJD(5) + t276 * t271 + (t310 * t271 - t321 * t274) * qJD(1), t304, t271 * t282 + (t281 * t267 - t240) * t274 + t277, t271 * t279 + (t281 + t315) * t306 + t277, t263, -t234 * r_i_i_C(1) + t235 * r_i_i_C(2); 0, 0, t275 - t299, t275, 0, t278 * r_i_i_C(2) + (-t289 + t297) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end