% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:52
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:52:18
	% EndTime: 2019-10-10 10:52:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:52:18
	% EndTime: 2019-10-10 10:52:18
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
	% StartTime: 2019-10-10 10:52:18
	% EndTime: 2019-10-10 10:52:19
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-10-10 10:52:18
	% EndTime: 2019-10-10 10:52:19
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (50->16), mult. (90->25), div. (0->0), fcn. (59->6), ass. (0->14)
	t27 = qJ(2) + pkin(11);
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
	% StartTime: 2019-10-10 10:52:18
	% EndTime: 2019-10-10 10:52:19
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (131->30), mult. (132->40), div. (0->0), fcn. (86->8), ass. (0->27)
	t49 = qJ(2) + pkin(11);
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
	% StartTime: 2019-10-10 10:52:19
	% EndTime: 2019-10-10 10:52:19
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (276->41), mult. (188->50), div. (0->0), fcn. (122->10), ass. (0->37)
	t58 = qJD(2) + qJD(4);
	t55 = qJD(5) + t58;
	t59 = qJ(2) + pkin(11);
	t56 = qJ(4) + t59;
	t52 = qJ(5) + t56;
	t49 = cos(t52);
	t83 = r_i_i_C(2) * t49;
	t48 = sin(t52);
	t85 = r_i_i_C(1) * t48;
	t68 = t83 + t85;
	t66 = t68 * t55;
	t50 = sin(t56);
	t86 = pkin(4) * t50;
	t88 = -t58 * t86 - t66;
	t70 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t59);
	t64 = t70 * qJD(2) + t88;
	t81 = t49 * t55;
	t75 = r_i_i_C(1) * t81;
	t51 = cos(t56);
	t80 = t51 * t58;
	t87 = -pkin(4) * t80 - t75;
	t84 = r_i_i_C(2) * t48;
	t82 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(3) + pkin(7);
	t61 = sin(qJ(1));
	t79 = qJD(1) * t61;
	t63 = cos(qJ(1));
	t78 = qJD(1) * t63;
	t74 = t55 * t84;
	t72 = qJD(1) * t83;
	t73 = t61 * t72 + t63 * t74 + t79 * t85;
	t69 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t59);
	t71 = t69 * qJD(2) + t87;
	t67 = -pkin(4) * t51 - r_i_i_C(1) * t49 - pkin(1) + t69 + t84;
	t65 = -t63 * t75 + t73;
	t44 = t61 * t74;
	t43 = t70 - t86;
	t1 = [t63 * qJD(3) - t64 * t61 + (-t82 * t61 + t67 * t63) * qJD(1), -t43 * t79 + t71 * t63 + t73, t78, (t50 * t79 - t63 * t80) * pkin(4) + t65, t65, 0; t61 * qJD(3) + t64 * t63 + (t67 * t61 + t82 * t63) * qJD(1), t44 + t71 * t61 + (t43 - t68) * t78, t79, t44 + t87 * t61 + (-t68 - t86) * t78, -t63 * t72 + t44 + (-t48 * t78 - t61 * t81) * r_i_i_C(1), 0; 0, t64, 0, t88, -t66, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:52:20
	% EndTime: 2019-10-10 10:52:20
	% DurationCPUTime: 0.48s
	% Computational Cost: add. (745->76), mult. (542->107), div. (0->0), fcn. (399->12), ass. (0->64)
	t278 = qJ(2) + pkin(11);
	t275 = qJ(4) + t278;
	t271 = qJ(5) + t275;
	t267 = sin(t271);
	t282 = cos(qJ(6));
	t338 = r_i_i_C(1) * t282 + pkin(5);
	t301 = t338 * t267;
	t279 = sin(qJ(6));
	t321 = qJD(6) * t282;
	t268 = cos(t271);
	t277 = qJD(2) + qJD(4);
	t274 = qJD(5) + t277;
	t329 = t268 * t274;
	t340 = t267 * t321 + t279 * t329;
	t335 = pkin(10) + r_i_i_C(3);
	t315 = t335 * t268;
	t298 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t278);
	t293 = t298 * qJD(2);
	t269 = sin(t275);
	t333 = pkin(4) * t277;
	t320 = t269 * t333;
	t332 = pkin(5) * t267;
	t339 = (t315 - t332) * t274 + t293 - t320;
	t322 = qJD(6) * t279;
	t308 = t267 * t322;
	t336 = r_i_i_C(1) * t308 + t340 * r_i_i_C(2);
	t270 = cos(t275);
	t300 = t338 * t268;
	t316 = t335 * t267;
	t287 = (-t300 - t316) * t274 - t270 * t333;
	t334 = pkin(4) * t269;
	t330 = r_i_i_C(2) * t279;
	t281 = sin(qJ(1));
	t328 = t274 * t281;
	t327 = t274 * t282;
	t284 = cos(qJ(1));
	t326 = t274 * t284;
	t325 = t282 * t284;
	t324 = qJD(1) * t281;
	t323 = qJD(1) * t284;
	t318 = t267 * t330;
	t317 = qJD(1) * t330;
	t314 = t335 * t281;
	t313 = t267 * t327;
	t303 = qJD(6) * t268 - qJD(1);
	t302 = qJD(1) * t268 - qJD(6);
	t299 = t338 * t284;
	t297 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t278);
	t296 = t336 * t284 + t324 * t301;
	t295 = t303 * t279;
	t294 = t284 * t267 * t317 + t336 * t281 + t323 * t315;
	t292 = -t315 - t318;
	t291 = -pkin(4) * t270 - pkin(5) * t268 - pkin(1) + t297 - t316;
	t290 = t267 * t326 + t302 * t281;
	t288 = t297 * qJD(2) + t287;
	t286 = -t268 * r_i_i_C(2) * t321 + (-t268 * t322 - t313) * r_i_i_C(1) + t335 * t329 + (-t332 + t318) * t274;
	t285 = t286 - t320;
	t276 = -pkin(9) - pkin(8) - qJ(3) - pkin(7);
	t252 = t298 - t334;
	t248 = -t302 * t325 + (t295 + t313) * t281;
	t247 = t303 * t282 * t281 + (-t267 * t328 + t302 * t284) * t279;
	t246 = t290 * t282 + t284 * t295;
	t245 = t290 * t279 - t303 * t325;
	t1 = [t248 * r_i_i_C(1) + t247 * r_i_i_C(2) + t284 * qJD(3) - t339 * t281 + (t276 * t281 + t291 * t284) * qJD(1), (-t252 + t292) * t324 + t288 * t284 + t296, t323, (t292 + t334) * t324 + t287 * t284 + t296, (-t281 * t317 - t335 * t326) * t267 + (-qJD(1) * t314 - t274 * t299) * t268 + t296, t245 * r_i_i_C(1) + t246 * r_i_i_C(2); -t246 * r_i_i_C(1) + t245 * r_i_i_C(2) + t281 * qJD(3) + t339 * t284 + (-t276 * t284 + t291 * t281) * qJD(1), (t252 - t301) * t323 + t288 * t281 + t294, t324, (-t301 - t334) * t323 + t287 * t281 + t294, -t300 * t328 + (-qJD(1) * t299 - t274 * t314) * t267 + t294, -t247 * r_i_i_C(1) + t248 * r_i_i_C(2); 0, t293 + t285, 0, t285, t286, (-t268 * t327 + t308) * r_i_i_C(2) - t340 * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end