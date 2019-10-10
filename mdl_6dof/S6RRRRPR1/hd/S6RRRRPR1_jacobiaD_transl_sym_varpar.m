% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:35
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:07
	% EndTime: 2019-10-10 12:35:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:07
	% EndTime: 2019-10-10 12:35:07
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
	% StartTime: 2019-10-10 12:35:07
	% EndTime: 2019-10-10 12:35:07
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 12:35:07
	% EndTime: 2019-10-10 12:35:07
	% DurationCPUTime: 0.14s
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
	t58 = r_i_i_C(3) + pkin(8) + pkin(7);
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
	t1 = [t63 * t41 + (-t58 * t41 + t48 * t43) * qJD(1), (t40 * t56 - t43 * t54) * pkin(2) + t46, t46, 0, 0, 0; -t63 * t43 + (t48 * t41 + t58 * t43) * qJD(1), t31 + (-pkin(2) * t54 - t53) * t41 + (-t49 - t62) * t55, -t43 * t51 + t31 + (-t36 * t55 - t41 * t57) * r_i_i_C(1), 0, 0, 0; 0, -t63, -t47, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:07
	% EndTime: 2019-10-10 12:35:07
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (202->36), mult. (170->47), div. (0->0), fcn. (109->8), ass. (0->37)
	t54 = sin(qJ(2));
	t51 = qJD(2) + qJD(3);
	t47 = qJD(4) + t51;
	t53 = qJ(2) + qJ(3);
	t50 = qJ(4) + t53;
	t46 = cos(t50);
	t77 = r_i_i_C(2) * t46;
	t45 = sin(t50);
	t79 = r_i_i_C(1) * t45;
	t62 = t77 + t79;
	t60 = t62 * t47;
	t48 = sin(t53);
	t80 = pkin(3) * t48;
	t58 = -t51 * t80 - t60;
	t73 = pkin(2) * qJD(2);
	t82 = -t54 * t73 + t58;
	t75 = t46 * t47;
	t68 = r_i_i_C(1) * t75;
	t49 = cos(t53);
	t74 = t49 * t51;
	t81 = -pkin(3) * t74 - t68;
	t78 = r_i_i_C(2) * t45;
	t76 = r_i_i_C(3) + pkin(9) + pkin(8) + pkin(7);
	t55 = sin(qJ(1));
	t72 = qJD(1) * t55;
	t57 = cos(qJ(1));
	t71 = qJD(1) * t57;
	t67 = t47 * t78;
	t65 = qJD(1) * t77;
	t66 = t55 * t65 + t57 * t67 + t72 * t79;
	t56 = cos(qJ(2));
	t63 = -t56 * t73 + t81;
	t61 = -pkin(2) * t56 - pkin(3) * t49 - r_i_i_C(1) * t46 - pkin(1) + t78;
	t59 = -t57 * t68 + t66;
	t44 = -pkin(2) * t54 - t80;
	t39 = t55 * t67;
	t1 = [-t82 * t55 + (-t76 * t55 + t61 * t57) * qJD(1), -t44 * t72 + t63 * t57 + t66, (t48 * t72 - t57 * t74) * pkin(3) + t59, t59, 0, 0; t82 * t57 + (t61 * t55 + t76 * t57) * qJD(1), t39 + t63 * t55 + (t44 - t62) * t71, t39 + t81 * t55 + (-t62 - t80) * t71, -t57 * t65 + t39 + (-t45 * t71 - t55 * t75) * r_i_i_C(1), 0, 0; 0, t82, t58, -t60, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:07
	% EndTime: 2019-10-10 12:35:07
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (323->43), mult. (210->51), div. (0->0), fcn. (136->10), ass. (0->40)
	t64 = sin(qJ(2));
	t63 = qJ(2) + qJ(3);
	t58 = sin(t63);
	t62 = qJD(2) + qJD(3);
	t57 = qJD(4) + t62;
	t61 = qJ(4) + t63;
	t54 = pkin(11) + t61;
	t52 = sin(t54);
	t53 = cos(t54);
	t72 = r_i_i_C(1) * t52 + r_i_i_C(2) * t53;
	t55 = sin(t61);
	t90 = pkin(4) * t55;
	t94 = t72 + t90;
	t69 = t94 * t57;
	t91 = pkin(3) * t62;
	t68 = -t58 * t91 - t69;
	t82 = pkin(2) * qJD(2);
	t95 = -t64 * t82 + t68;
	t59 = cos(t63);
	t87 = r_i_i_C(1) * t53;
	t78 = t57 * t87;
	t56 = cos(t61);
	t83 = t56 * t57;
	t73 = -pkin(4) * t83 - t59 * t91 - t78;
	t93 = -pkin(4) * t56 - t87;
	t86 = r_i_i_C(2) * t52;
	t84 = r_i_i_C(3) + qJ(5) + pkin(9) + pkin(8) + pkin(7);
	t65 = sin(qJ(1));
	t81 = qJD(1) * t65;
	t67 = cos(qJ(1));
	t80 = qJD(1) * t67;
	t77 = t57 * t86;
	t76 = t67 * t77 + t72 * t81;
	t66 = cos(qJ(2));
	t74 = -t66 * t82 + t73;
	t49 = -pkin(3) * t58 - t90;
	t71 = -t66 * pkin(2) - pkin(3) * t59 - pkin(1) + t86 + t93;
	t47 = t65 * t77;
	t46 = -t64 * pkin(2) + t49;
	t1 = [t67 * qJD(5) - t95 * t65 + (-t84 * t65 + t71 * t67) * qJD(1), -t46 * t81 + t74 * t67 + t76, -t49 * t81 + t73 * t67 + t76, -t67 * t78 + (t55 * t81 - t67 * t83) * pkin(4) + t76, t80, 0; t65 * qJD(5) + t95 * t67 + (t71 * t65 + t84 * t67) * qJD(1), t47 + t74 * t65 + (t46 - t72) * t80, t47 + t73 * t65 + (t49 - t72) * t80, t93 * t65 * t57 - t80 * t94 + t47, t81, 0; 0, t95, t68, -t69, 0, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:08
	% EndTime: 2019-10-10 12:35:09
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (792->76), mult. (564->97), div. (0->0), fcn. (413->12), ass. (0->65)
	t283 = qJ(2) + qJ(3);
	t281 = qJ(4) + t283;
	t274 = pkin(11) + t281;
	t272 = sin(t274);
	t287 = cos(qJ(6));
	t339 = r_i_i_C(1) * t287 + pkin(5);
	t304 = t339 * t272;
	t284 = sin(qJ(6));
	t322 = qJD(6) * t272;
	t273 = cos(t274);
	t282 = qJD(2) + qJD(3);
	t277 = qJD(4) + t282;
	t327 = t273 * t277;
	t342 = t284 * t327 + t287 * t322;
	t336 = pkin(10) + r_i_i_C(3);
	t316 = t336 * t273;
	t319 = r_i_i_C(2) * t272 * t284;
	t341 = t316 + t319;
	t275 = sin(t281);
	t285 = sin(qJ(2));
	t329 = pkin(2) * qJD(2);
	t318 = t285 * t329;
	t278 = sin(t283);
	t335 = pkin(3) * t282;
	t320 = t278 * t335;
	t332 = pkin(4) * t277;
	t340 = -t275 * t332 + (-pkin(5) * t272 + t316) * t277 - t318 - t320;
	t311 = t284 * t322;
	t337 = r_i_i_C(1) * t311 + t342 * r_i_i_C(2);
	t276 = cos(t281);
	t279 = cos(t283);
	t317 = t336 * t272;
	t297 = -t273 * t339 - t317;
	t293 = -t276 * t332 + t297 * t277 - t279 * t335;
	t334 = pkin(4) * t275;
	t333 = pkin(4) * t276;
	t328 = t272 * t277;
	t326 = t277 * t287;
	t289 = cos(qJ(1));
	t325 = t287 * t289;
	t286 = sin(qJ(1));
	t324 = qJD(1) * t286;
	t323 = qJD(1) * t289;
	t321 = qJD(6) * t273;
	t306 = -qJD(1) + t321;
	t305 = qJD(1) * t273 - qJD(6);
	t266 = -pkin(3) * t278 - t334;
	t303 = t337 * t289 + t324 * t304;
	t302 = t306 * t284;
	t301 = t337 * t286 + t341 * t323;
	t288 = cos(qJ(2));
	t299 = -t288 * pkin(2) - pkin(3) * t279 - pkin(5) * t273 - pkin(1) - t317 - t333;
	t298 = -t304 - t334;
	t296 = t305 * t286 + t289 * t328;
	t294 = -t288 * t329 + t293;
	t292 = t277 * (t297 - t333);
	t291 = (-r_i_i_C(1) * t284 - r_i_i_C(2) * t287) * t321 + t336 * t327 + (t319 + t298) * t277;
	t290 = t291 - t320;
	t280 = -qJ(5) - pkin(9) - pkin(8) - pkin(7);
	t258 = -t285 * pkin(2) + t266;
	t251 = -t305 * t325 + (t272 * t326 + t302) * t286;
	t250 = t306 * t287 * t286 + (-t286 * t328 + t305 * t289) * t284;
	t249 = t296 * t287 + t289 * t302;
	t248 = t296 * t284 - t306 * t325;
	t1 = [t251 * r_i_i_C(1) + t250 * r_i_i_C(2) + t289 * qJD(5) - t340 * t286 + (t280 * t286 + t299 * t289) * qJD(1), (-t258 - t341) * t324 + t294 * t289 + t303, (-t266 - t341) * t324 + t293 * t289 + t303, (-t341 + t334) * t324 + t289 * t292 + t303, t323, t248 * r_i_i_C(1) + t249 * r_i_i_C(2); -t249 * r_i_i_C(1) + t248 * r_i_i_C(2) + t286 * qJD(5) + t340 * t289 + (-t280 * t289 + t299 * t286) * qJD(1), (t258 - t304) * t323 + t294 * t286 + t301, (t266 - t304) * t323 + t293 * t286 + t301, t286 * t292 + t298 * t323 + t301, t324, -t250 * r_i_i_C(1) + t251 * r_i_i_C(2); 0, t290 - t318, t290, t291, 0, (-t273 * t326 + t311) * r_i_i_C(2) - t342 * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end