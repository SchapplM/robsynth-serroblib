% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:55
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:55:18
	% EndTime: 2019-10-10 11:55:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:55:18
	% EndTime: 2019-10-10 11:55:18
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
	% StartTime: 2019-10-10 11:55:18
	% EndTime: 2019-10-10 11:55:18
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
	% StartTime: 2019-10-10 11:55:18
	% EndTime: 2019-10-10 11:55:18
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
	% StartTime: 2019-10-10 11:55:18
	% EndTime: 2019-10-10 11:55:18
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (149->33), mult. (146->42), div. (0->0), fcn. (95->8), ass. (0->32)
	t52 = sin(qJ(2));
	t50 = qJD(2) + qJD(3);
	t51 = qJ(2) + qJ(3);
	t46 = pkin(11) + t51;
	t44 = sin(t46);
	t45 = cos(t46);
	t59 = r_i_i_C(1) * t44 + r_i_i_C(2) * t45;
	t47 = sin(t51);
	t75 = pkin(3) * t47;
	t78 = t59 + t75;
	t56 = t78 * t50;
	t67 = pkin(2) * qJD(2);
	t76 = t52 * t67 + t56;
	t48 = cos(t51);
	t72 = r_i_i_C(1) * t45;
	t77 = -pkin(3) * t48 - t72;
	t71 = r_i_i_C(2) * t44;
	t69 = r_i_i_C(3) + qJ(4) + pkin(8) + pkin(7);
	t68 = t48 * t50;
	t53 = sin(qJ(1));
	t66 = qJD(1) * t53;
	t55 = cos(qJ(1));
	t65 = qJD(1) * t55;
	t64 = t50 * t72;
	t63 = t50 * t71;
	t62 = t55 * t63 + t59 * t66;
	t54 = cos(qJ(2));
	t60 = -pkin(3) * t68 - t54 * t67 - t64;
	t58 = -t54 * pkin(2) - pkin(1) + t71 + t77;
	t43 = -t52 * pkin(2) - t75;
	t38 = t53 * t63;
	t1 = [t55 * qJD(4) + t76 * t53 + (-t69 * t53 + t58 * t55) * qJD(1), -t43 * t66 + t60 * t55 + t62, -t55 * t64 + (t47 * t66 - t55 * t68) * pkin(3) + t62, t65, 0, 0; t53 * qJD(4) - t76 * t55 + (t58 * t53 + t69 * t55) * qJD(1), t38 + t60 * t53 + (t43 - t59) * t65, t50 * t53 * t77 - t65 * t78 + t38, t66, 0, 0; 0, -t76, -t56, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:55:18
	% EndTime: 2019-10-10 11:55:18
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (294->41), mult. (196->49), div. (0->0), fcn. (127->10), ass. (0->35)
	t63 = sin(qJ(2));
	t62 = qJ(2) + qJ(3);
	t57 = pkin(11) + t62;
	t48 = -pkin(3) * sin(t62) - pkin(4) * sin(t57);
	t61 = qJD(2) + qJD(3);
	t56 = qJD(5) + t61;
	t55 = qJ(5) + t57;
	t52 = cos(t55);
	t85 = r_i_i_C(2) * t52;
	t51 = sin(t55);
	t87 = r_i_i_C(1) * t51;
	t71 = t85 + t87;
	t68 = t71 * t56;
	t67 = t48 * t61 - t68;
	t82 = pkin(2) * qJD(2);
	t88 = -t63 * t82 + t67;
	t72 = -pkin(3) * cos(t62) - pkin(4) * cos(t57);
	t83 = t52 * t56;
	t79 = r_i_i_C(1) * t83;
	t73 = t72 * t61 - t79;
	t86 = r_i_i_C(2) * t51;
	t84 = r_i_i_C(3) + pkin(9) + qJ(4) + pkin(8) + pkin(7);
	t64 = sin(qJ(1));
	t81 = qJD(1) * t64;
	t66 = cos(qJ(1));
	t80 = qJD(1) * t66;
	t78 = t56 * t86;
	t76 = qJD(1) * t85;
	t77 = t64 * t76 + t66 * t78 + t81 * t87;
	t65 = cos(qJ(2));
	t74 = -t65 * t82 + t73;
	t70 = -t65 * pkin(2) - r_i_i_C(1) * t52 - pkin(1) + t72 + t86;
	t46 = t64 * t78;
	t45 = -t63 * pkin(2) + t48;
	t1 = [t66 * qJD(4) - t88 * t64 + (-t84 * t64 + t70 * t66) * qJD(1), -t45 * t81 + t74 * t66 + t77, -t48 * t81 + t73 * t66 + t77, t80, -t66 * t79 + t77, 0; t64 * qJD(4) + t88 * t66 + (t70 * t64 + t84 * t66) * qJD(1), t46 + t74 * t64 + (t45 - t71) * t80, t46 + t73 * t64 + (t48 - t71) * t80, t81, -t66 * t76 + t46 + (-t51 * t80 - t64 * t83) * r_i_i_C(1), 0; 0, t88, t67, 0, -t68, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:55:19
	% EndTime: 2019-10-10 11:55:20
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (763->76), mult. (550->107), div. (0->0), fcn. (404->12), ass. (0->63)
	t281 = qJ(2) + qJ(3);
	t276 = pkin(11) + t281;
	t274 = qJ(5) + t276;
	t270 = sin(t274);
	t285 = cos(qJ(6));
	t338 = r_i_i_C(1) * t285 + pkin(5);
	t303 = t338 * t270;
	t282 = sin(qJ(6));
	t322 = qJD(6) * t285;
	t271 = cos(t274);
	t280 = qJD(2) + qJD(3);
	t275 = qJD(5) + t280;
	t330 = t271 * t275;
	t340 = t270 * t322 + t282 * t330;
	t335 = pkin(10) + r_i_i_C(3);
	t317 = t335 * t271;
	t264 = -pkin(3) * sin(t281) - pkin(4) * sin(t276);
	t296 = t264 * t280;
	t283 = sin(qJ(2));
	t331 = pkin(2) * qJD(2);
	t319 = t283 * t331;
	t334 = pkin(5) * t270;
	t339 = (t317 - t334) * t275 + t296 - t319;
	t323 = qJD(6) * t282;
	t310 = t270 * t323;
	t336 = r_i_i_C(1) * t310 + r_i_i_C(2) * t340;
	t300 = -pkin(3) * cos(t281) - pkin(4) * cos(t276);
	t302 = t338 * t271;
	t318 = t335 * t270;
	t290 = t300 * t280 + (-t302 - t318) * t275;
	t332 = r_i_i_C(2) * t282;
	t284 = sin(qJ(1));
	t329 = t275 * t284;
	t328 = t275 * t285;
	t287 = cos(qJ(1));
	t327 = t275 * t287;
	t326 = t285 * t287;
	t325 = qJD(1) * t284;
	t324 = qJD(1) * t287;
	t321 = t270 * t332;
	t320 = qJD(1) * t332;
	t316 = t335 * t284;
	t315 = t270 * t328;
	t305 = qJD(6) * t271 - qJD(1);
	t304 = qJD(1) * t271 - qJD(6);
	t301 = t338 * t287;
	t299 = t287 * t336 + t303 * t325;
	t298 = t305 * t282;
	t297 = t287 * t270 * t320 + t284 * t336 + t317 * t324;
	t295 = -t317 - t321;
	t286 = cos(qJ(2));
	t294 = -pkin(2) * t286 - pkin(5) * t271 - pkin(1) + t300 - t318;
	t293 = t270 * t327 + t284 * t304;
	t291 = -t286 * t331 + t290;
	t289 = -t271 * r_i_i_C(2) * t322 + (-t271 * t323 - t315) * r_i_i_C(1) + t335 * t330 + (-t334 + t321) * t275;
	t288 = t296 + t289;
	t279 = -pkin(9) - qJ(4) - pkin(8) - pkin(7);
	t256 = -pkin(2) * t283 + t264;
	t249 = -t304 * t326 + (t298 + t315) * t284;
	t248 = t305 * t285 * t284 + (-t270 * t329 + t287 * t304) * t282;
	t247 = t285 * t293 + t287 * t298;
	t246 = t282 * t293 - t305 * t326;
	t1 = [t249 * r_i_i_C(1) + t248 * r_i_i_C(2) + t287 * qJD(4) - t339 * t284 + (t279 * t284 + t287 * t294) * qJD(1), (-t256 + t295) * t325 + t291 * t287 + t299, (-t264 + t295) * t325 + t290 * t287 + t299, t324, (-t284 * t320 - t327 * t335) * t270 + (-qJD(1) * t316 - t275 * t301) * t271 + t299, r_i_i_C(1) * t246 + r_i_i_C(2) * t247; -t247 * r_i_i_C(1) + t246 * r_i_i_C(2) + t284 * qJD(4) + t339 * t287 + (-t279 * t287 + t284 * t294) * qJD(1), (t256 - t303) * t324 + t291 * t284 + t297, (t264 - t303) * t324 + t290 * t284 + t297, t325, -t302 * t329 + (-qJD(1) * t301 - t275 * t316) * t270 + t297, -r_i_i_C(1) * t248 + r_i_i_C(2) * t249; 0, t288 - t319, t288, 0, t289, (-t271 * t328 + t310) * r_i_i_C(2) - t340 * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end