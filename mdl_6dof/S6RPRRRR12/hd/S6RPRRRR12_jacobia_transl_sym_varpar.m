% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRRR12
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:16
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR12_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR12_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR12_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobia_transl_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (13->12), mult. (28->19), div. (0->0), fcn. (35->6), ass. (0->12)
	t1 = sin(pkin(14));
	t5 = sin(qJ(1));
	t11 = t5 * t1;
	t3 = cos(pkin(14));
	t10 = t5 * t3;
	t6 = cos(qJ(1));
	t9 = t6 * t1;
	t8 = t6 * t3;
	t2 = sin(pkin(6));
	t7 = t2 * (r_i_i_C(3) + qJ(2));
	t4 = cos(pkin(6));
	t12 = [(-t4 * t9 - t10) * r_i_i_C(1) + (-t4 * t8 + t11) * r_i_i_C(2) - t5 * pkin(1) + t6 * t7, t5 * t2, 0, 0, 0, 0; (-t11 * t4 + t8) * r_i_i_C(1) + (-t10 * t4 - t9) * r_i_i_C(2) + t6 * pkin(1) + t5 * t7, -t6 * t2, 0, 0, 0, 0; 0, t4, 0, 0, 0, 0;];
	Ja_transl = t12;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (62->32), mult. (166->58), div. (0->0), fcn. (213->10), ass. (0->29)
	t34 = r_i_i_C(3) + pkin(10);
	t13 = cos(pkin(7));
	t15 = sin(qJ(3));
	t17 = cos(qJ(3));
	t10 = sin(pkin(7));
	t11 = sin(pkin(6));
	t18 = cos(qJ(1));
	t26 = t18 * t11;
	t24 = t10 * t26;
	t12 = cos(pkin(14));
	t14 = cos(pkin(6));
	t29 = t14 * t18;
	t16 = sin(qJ(1));
	t9 = sin(pkin(14));
	t32 = t16 * t9;
	t3 = -t12 * t29 + t32;
	t27 = t16 * t12;
	t4 = t9 * t29 + t27;
	t33 = (t13 * t3 + t24) * t17 + t4 * t15;
	t30 = t13 * t15;
	t28 = t16 * t11;
	t25 = t11 * qJ(2);
	t5 = -t14 * t27 - t18 * t9;
	t22 = t10 * t28 + t13 * t5;
	t19 = t15 * t24 - t17 * t4 + t3 * t30;
	t6 = t12 * t18 - t14 * t32;
	t2 = t22 * t15 + t17 * t6;
	t1 = -t15 * t6 + t22 * t17;
	t7 = [t19 * r_i_i_C(1) + t33 * r_i_i_C(2) - t4 * pkin(2) - t16 * pkin(1) + t18 * t25 + t34 * (-t3 * t10 + t13 * t26), t28, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; t18 * pkin(1) + t6 * pkin(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t16 * t25 + t34 * (-t10 * t5 + t13 * t28), -t26, -t33 * r_i_i_C(1) + t19 * r_i_i_C(2), 0, 0, 0; 0, t14, (t17 * r_i_i_C(1) - t15 * r_i_i_C(2)) * t14 * t10 + ((t12 * t13 * t17 - t15 * t9) * r_i_i_C(1) + (-t12 * t30 - t17 * t9) * r_i_i_C(2)) * t11, 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (260->63), mult. (734->112), div. (0->0), fcn. (961->14), ass. (0->49)
	t29 = sin(qJ(4));
	t32 = cos(qJ(4));
	t28 = cos(pkin(6));
	t25 = cos(pkin(14));
	t34 = cos(qJ(1));
	t46 = t34 * t25;
	t21 = sin(pkin(14));
	t31 = sin(qJ(1));
	t51 = t31 * t21;
	t15 = -t28 * t46 + t51;
	t23 = sin(pkin(7));
	t27 = cos(pkin(7));
	t24 = sin(pkin(6));
	t47 = t34 * t24;
	t11 = -t15 * t23 + t27 * t47;
	t22 = sin(pkin(8));
	t26 = cos(pkin(8));
	t48 = t34 * t21;
	t49 = t31 * t25;
	t16 = t28 * t48 + t49;
	t30 = sin(qJ(3));
	t33 = cos(qJ(3));
	t44 = t23 * t47;
	t5 = (t15 * t27 + t44) * t33 + t16 * t30;
	t41 = t11 * t22 + t26 * t5;
	t52 = t27 * t30;
	t6 = t15 * t52 - t16 * t33 + t30 * t44;
	t62 = t41 * t29 + t6 * t32;
	t61 = -t6 * t29 + t41 * t32;
	t58 = pkin(11) + r_i_i_C(3);
	t55 = t23 * t28;
	t54 = t26 * t29;
	t53 = t26 * t32;
	t50 = t31 * t24;
	t45 = t24 * qJ(2);
	t43 = t58 * t22;
	t17 = -t28 * t49 - t48;
	t13 = -t17 * t23 + t27 * t50;
	t18 = -t28 * t51 + t46;
	t35 = t17 * t27 + t23 * t50;
	t7 = -t18 * t30 + t35 * t33;
	t38 = t13 * t22 + t26 * t7;
	t9 = t33 * t55 + (t25 * t27 * t33 - t21 * t30) * t24;
	t37 = (-t24 * t25 * t23 + t28 * t27) * t22 + t26 * t9;
	t10 = t30 * t55 + (t21 * t33 + t25 * t52) * t24;
	t8 = t18 * t33 + t35 * t30;
	t2 = t38 * t29 + t8 * t32;
	t1 = -t8 * t29 + t38 * t32;
	t3 = [t62 * r_i_i_C(1) + t61 * r_i_i_C(2) + t6 * pkin(3) - t16 * pkin(2) - t31 * pkin(1) + t34 * t45 + t11 * pkin(10) + t58 * (t11 * t26 - t5 * t22), t50, (t7 * t32 - t8 * t54) * r_i_i_C(1) + (-t7 * t29 - t8 * t53) * r_i_i_C(2) + t7 * pkin(3) + t8 * t43, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t34 * pkin(1) + t18 * pkin(2) + t8 * pkin(3) + t13 * pkin(10) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t31 * t45 + t58 * (t13 * t26 - t7 * t22), -t47, (-t32 * t5 + t54 * t6) * r_i_i_C(1) + (t29 * t5 + t53 * t6) * r_i_i_C(2) - t5 * pkin(3) - t6 * t43, -t61 * r_i_i_C(1) + t62 * r_i_i_C(2), 0, 0; 0, t28, (t32 * r_i_i_C(1) - t29 * r_i_i_C(2) + pkin(3)) * t9 + ((-t29 * r_i_i_C(1) - t32 * r_i_i_C(2)) * t26 + t43) * t10, (-t10 * t29 + t37 * t32) * r_i_i_C(1) + (-t10 * t32 - t37 * t29) * r_i_i_C(2), 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.76s
	% Computational Cost: add. (654->92), mult. (1858->160), div. (0->0), fcn. (2457->16), ass. (0->66)
	t47 = cos(pkin(6));
	t44 = cos(pkin(14));
	t55 = cos(qJ(1));
	t66 = t55 * t44;
	t40 = sin(pkin(14));
	t51 = sin(qJ(1));
	t71 = t51 * t40;
	t34 = -t47 * t66 + t71;
	t68 = t55 * t40;
	t69 = t51 * t44;
	t35 = t47 * t68 + t69;
	t46 = cos(pkin(7));
	t50 = sin(qJ(3));
	t54 = cos(qJ(3));
	t42 = sin(pkin(7));
	t43 = sin(pkin(6));
	t67 = t55 * t43;
	t64 = t42 * t67;
	t24 = (t34 * t46 + t64) * t54 + t35 * t50;
	t31 = -t34 * t42 + t46 * t67;
	t41 = sin(pkin(8));
	t45 = cos(pkin(8));
	t15 = t24 * t41 - t31 * t45;
	t48 = sin(qJ(5));
	t52 = cos(qJ(5));
	t72 = t46 * t50;
	t25 = t34 * t72 - t35 * t54 + t50 * t64;
	t49 = sin(qJ(4));
	t53 = cos(qJ(4));
	t62 = t24 * t45 + t31 * t41;
	t6 = t25 * t53 + t62 * t49;
	t89 = t15 * t52 + t6 * t48;
	t88 = -t15 * t48 + t6 * t52;
	t85 = t25 * t49 - t62 * t53;
	t80 = r_i_i_C(3) + pkin(12);
	t79 = t41 * pkin(11);
	t77 = t41 * t48;
	t76 = t41 * t52;
	t75 = t42 * t47;
	t74 = t45 * t49;
	t73 = t45 * t53;
	t70 = t51 * t43;
	t65 = t43 * qJ(2);
	t29 = t54 * t75 + (t44 * t46 * t54 - t40 * t50) * t43;
	t33 = -t43 * t44 * t42 + t47 * t46;
	t61 = t29 * t45 + t33 * t41;
	t60 = t52 * r_i_i_C(1) - t48 * r_i_i_C(2) + pkin(4);
	t36 = -t47 * t69 - t68;
	t58 = -t36 * t42 + t46 * t70;
	t57 = t36 * t46 + t42 * t70;
	t56 = t58 * t41;
	t37 = -t47 * t71 + t66;
	t26 = -t37 * t50 + t57 * t54;
	t17 = -t26 * t41 + t58 * t45;
	t30 = t50 * t75 + (t40 * t54 + t44 * t72) * t43;
	t27 = t37 * t54 + t57 * t50;
	t21 = -t29 * t41 + t33 * t45;
	t19 = t29 * t53 - t30 * t74;
	t14 = t30 * t53 + t61 * t49;
	t12 = t26 * t53 - t27 * t74;
	t10 = -t24 * t53 + t25 * t74;
	t8 = t27 * t53 + (t26 * t45 + t56) * t49;
	t7 = -t26 * t73 + t27 * t49 - t53 * t56;
	t2 = t17 * t48 + t8 * t52;
	t1 = t17 * t52 - t8 * t48;
	t3 = [-t51 * pkin(1) - t35 * pkin(2) + t25 * pkin(3) + t6 * pkin(4) + t31 * pkin(10) - t15 * pkin(11) + t88 * r_i_i_C(1) - t89 * r_i_i_C(2) + t55 * t65 + t80 * t85, t70, (t12 * t52 + t27 * t77) * r_i_i_C(1) + (-t12 * t48 + t27 * t76) * r_i_i_C(2) + t12 * pkin(4) + t26 * pkin(3) + t27 * t79 + t80 * (t26 * t49 + t27 * t73), -t60 * t7 + t80 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t55 * pkin(1) + t37 * pkin(2) + t27 * pkin(3) + t8 * pkin(4) + t58 * pkin(10) + t17 * pkin(11) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t51 * t65 + t80 * t7, -t67, (t10 * t52 - t25 * t77) * r_i_i_C(1) + (-t10 * t48 - t25 * t76) * r_i_i_C(2) + t10 * pkin(4) - t24 * pkin(3) - t25 * t79 + t80 * (-t24 * t49 - t25 * t73), -t6 * t80 + t60 * t85, t89 * r_i_i_C(1) + t88 * r_i_i_C(2), 0; 0, t47, (t19 * t52 + t30 * t77) * r_i_i_C(1) + (-t19 * t48 + t30 * t76) * r_i_i_C(2) + t19 * pkin(4) + t29 * pkin(3) + t30 * t79 + t80 * (t29 * t49 + t30 * t73), t80 * t14 + t60 * (-t30 * t49 + t61 * t53), (-t14 * t48 + t21 * t52) * r_i_i_C(1) + (-t14 * t52 - t21 * t48) * r_i_i_C(2), 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:47
	% DurationCPUTime: 1.39s
	% Computational Cost: add. (1425->125), mult. (4066->209), div. (0->0), fcn. (5408->18), ass. (0->86)
	t107 = cos(qJ(4));
	t106 = sin(qJ(1));
	t50 = sin(pkin(14));
	t61 = cos(qJ(1));
	t97 = cos(pkin(14));
	t99 = cos(pkin(6));
	t88 = t99 * t97;
	t47 = t106 * t50 - t61 * t88;
	t52 = sin(pkin(6));
	t96 = sin(pkin(7));
	t91 = t52 * t96;
	t98 = cos(pkin(7));
	t113 = t98 * t47 + t61 * t91;
	t93 = t50 * t99;
	t48 = t106 * t97 + t61 * t93;
	t57 = sin(qJ(3));
	t60 = cos(qJ(3));
	t110 = t113 * t60 + t48 * t57;
	t51 = sin(pkin(8));
	t92 = t52 * t98;
	t82 = t47 * t96 - t61 * t92;
	t114 = t82 * t51;
	t37 = -t113 * t57 + t48 * t60;
	t53 = cos(pkin(8));
	t56 = sin(qJ(4));
	t16 = t37 * t107 + (-t110 * t53 + t114) * t56;
	t30 = t110 * t51 + t82 * t53;
	t55 = sin(qJ(5));
	t59 = cos(qJ(5));
	t4 = t16 * t59 + t30 * t55;
	t54 = sin(qJ(6));
	t122 = t4 * t54;
	t58 = cos(qJ(6));
	t121 = t4 * t58;
	t120 = -t16 * t55 + t30 * t59;
	t116 = t107 * t114 - t37 * t56;
	t109 = r_i_i_C(3) + pkin(13);
	t80 = -t106 * t88 - t61 * t50;
	t70 = t106 * t91 + t80 * t98;
	t79 = -t106 * t93 + t61 * t97;
	t65 = t79 * t57 - t70 * t60;
	t69 = t106 * t92 - t80 * t96;
	t62 = t65 * t51 + t69 * t53;
	t108 = pkin(11) * t51;
	t103 = t51 * t55;
	t102 = t51 * t59;
	t101 = t53 * t56;
	t100 = t61 * t52;
	t95 = t53 * t107;
	t94 = t106 * t52;
	t87 = t99 * t96;
	t86 = t98 * t97;
	t85 = t58 * r_i_i_C(1) - t54 * r_i_i_C(2) + pkin(5);
	t84 = t54 * r_i_i_C(1) + t58 * r_i_i_C(2) + pkin(12);
	t77 = -t109 * t55 - t85 * t59 - pkin(4);
	t76 = -t97 * t91 + t99 * t98;
	t73 = t76 * t51;
	t72 = t110 * t107;
	t71 = t60 * t87 + (-t50 * t57 + t60 * t86) * t52;
	t68 = t71 * t107;
	t67 = t69 * t51;
	t63 = t65 * t107;
	t44 = t57 * t87 + (t50 * t60 + t57 * t86) * t52;
	t40 = t70 * t57 + t79 * t60;
	t36 = -t71 * t51 + t76 * t53;
	t33 = -t44 * t101 + t68;
	t32 = t44 * t95 + t71 * t56;
	t28 = t44 * t107 + (t71 * t53 + t73) * t56;
	t27 = -t107 * t73 + t44 * t56 - t53 * t68;
	t26 = t44 * t103 + t33 * t59;
	t24 = -t40 * t101 - t63;
	t23 = t40 * t95 - t65 * t56;
	t22 = -t37 * t101 - t72;
	t21 = -t110 * t56 + t37 * t95;
	t20 = t40 * t107 + (-t65 * t53 + t67) * t56;
	t19 = -t107 * t67 + t40 * t56 + t53 * t63;
	t17 = -t110 * t95 + t116;
	t15 = t53 * t72 - t116;
	t14 = t28 * t59 + t36 * t55;
	t12 = t40 * t103 + t24 * t59;
	t10 = t37 * t103 + t22 * t59;
	t8 = t20 * t59 + t55 * t62;
	t7 = t20 * t55 - t59 * t62;
	t2 = t19 * t54 + t8 * t58;
	t1 = t19 * t58 - t8 * t54;
	t3 = [(t17 * t54 - t121) * r_i_i_C(1) + (t17 * t58 + t122) * r_i_i_C(2) - t4 * pkin(5) - t16 * pkin(4) + t17 * pkin(12) - t37 * pkin(3) - t48 * pkin(2) - t106 * pkin(1) + qJ(2) * t100 + t109 * t120 - t30 * pkin(11) - t82 * pkin(10), t94, (t12 * t58 + t23 * t54) * r_i_i_C(1) + (-t12 * t54 + t23 * t58) * r_i_i_C(2) + t12 * pkin(5) + t24 * pkin(4) + t23 * pkin(12) - t65 * pkin(3) + t40 * t108 + t109 * (-t40 * t102 + t24 * t55), t77 * t19 + t84 * t20, t109 * t8 - t85 * t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t61 * pkin(1) + t79 * pkin(2) + t40 * pkin(3) + t20 * pkin(4) + t8 * pkin(5) + t69 * pkin(10) + t62 * pkin(11) + t19 * pkin(12) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + qJ(2) * t94 + t109 * t7, -t100, (t10 * t58 + t21 * t54) * r_i_i_C(1) + (-t10 * t54 + t21 * t58) * r_i_i_C(2) + t10 * pkin(5) + t22 * pkin(4) + t21 * pkin(12) - t110 * pkin(3) + t37 * t108 + t109 * (-t37 * t102 + t22 * t55), t77 * t15 + t84 * t16, t109 * t4 + t85 * t120, (t15 * t58 - t122) * r_i_i_C(1) + (-t15 * t54 - t121) * r_i_i_C(2); 0, t99, (t26 * t58 + t32 * t54) * r_i_i_C(1) + (-t26 * t54 + t32 * t58) * r_i_i_C(2) + t26 * pkin(5) + t33 * pkin(4) + t32 * pkin(12) + t71 * pkin(3) + t44 * t108 + t109 * (-t44 * t102 + t33 * t55), t77 * t27 + t84 * t28, t109 * t14 + t85 * (-t28 * t55 + t36 * t59), (-t14 * t54 + t27 * t58) * r_i_i_C(1) + (-t14 * t58 - t27 * t54) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end