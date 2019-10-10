% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:35
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR10_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR10_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobia_transl_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (19->14), mult. (46->24), div. (0->0), fcn. (54->6), ass. (0->16)
	t7 = sin(qJ(2));
	t8 = sin(qJ(1));
	t15 = t8 * t7;
	t9 = cos(qJ(2));
	t14 = t8 * t9;
	t10 = cos(qJ(1));
	t13 = t10 * t7;
	t12 = t10 * t9;
	t5 = sin(pkin(6));
	t11 = (pkin(10) + r_i_i_C(3)) * t5;
	t6 = cos(pkin(6));
	t4 = -t6 * t15 + t12;
	t3 = -t6 * t14 - t13;
	t2 = -t6 * t13 - t14;
	t1 = -t6 * t12 + t15;
	t16 = [-t8 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t10 * t11, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0, 0, 0; t10 * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t11 * t8, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0, 0, 0; 0, (r_i_i_C(1) * t9 - r_i_i_C(2) * t7) * t5, 0, 0, 0, 0;];
	Ja_transl = t16;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (93->46), mult. (252->84), div. (0->0), fcn. (316->10), ass. (0->33)
	t39 = pkin(11) + r_i_i_C(3);
	t12 = cos(pkin(7));
	t13 = sin(qJ(3));
	t16 = cos(qJ(3));
	t10 = sin(pkin(7));
	t11 = sin(pkin(6));
	t18 = cos(qJ(1));
	t35 = t11 * t18;
	t27 = t10 * t35;
	t14 = sin(qJ(2));
	t15 = sin(qJ(1));
	t17 = cos(qJ(2));
	t28 = cos(pkin(6));
	t24 = t18 * t28;
	t3 = t15 * t14 - t17 * t24;
	t4 = t14 * t24 + t15 * t17;
	t38 = (t12 * t3 + t27) * t16 + t4 * t13;
	t36 = t11 * t15;
	t34 = t12 * t13;
	t33 = t12 * t16;
	t32 = t13 * t14;
	t31 = t13 * t17;
	t30 = t14 * t16;
	t29 = t16 * t17;
	t26 = t10 * t39;
	t25 = t15 * t28;
	t5 = -t18 * t14 - t17 * t25;
	t22 = t10 * t36 + t12 * t5;
	t19 = t13 * t27 - t4 * t16 + t3 * t34;
	t6 = -t14 * t25 + t18 * t17;
	t2 = t22 * t13 + t6 * t16;
	t1 = -t6 * t13 + t22 * t16;
	t7 = [t19 * r_i_i_C(1) + t38 * r_i_i_C(2) - t4 * pkin(2) - t15 * pkin(1) + pkin(10) * t35 + t39 * (-t3 * t10 + t12 * t35), (t5 * t16 - t6 * t34) * r_i_i_C(1) + (-t5 * t13 - t6 * t33) * r_i_i_C(2) + t5 * pkin(2) + t6 * t26, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; t18 * pkin(1) + t6 * pkin(2) + pkin(10) * t36 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t39 * (-t5 * t10 + t12 * t36), (-t3 * t16 - t4 * t34) * r_i_i_C(1) + (t3 * t13 - t4 * t33) * r_i_i_C(2) - t3 * pkin(2) + t4 * t26, -t38 * r_i_i_C(1) + t19 * r_i_i_C(2), 0, 0, 0; 0, ((-t12 * t32 + t29) * r_i_i_C(1) + (-t12 * t30 - t31) * r_i_i_C(2) + t17 * pkin(2) + t14 * t26) * t11, ((t12 * t29 - t32) * r_i_i_C(1) + (-t12 * t31 - t30) * r_i_i_C(2)) * t11 + (t16 * r_i_i_C(1) - t13 * r_i_i_C(2)) * t10 * t28, 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.71s
	% Computational Cost: add. (356->91), mult. (1016->164), div. (0->0), fcn. (1316->14), ass. (0->63)
	t33 = sin(qJ(4));
	t37 = cos(qJ(4));
	t35 = sin(qJ(2));
	t36 = sin(qJ(1));
	t39 = cos(qJ(2));
	t40 = cos(qJ(1));
	t60 = cos(pkin(6));
	t56 = t40 * t60;
	t21 = t36 * t35 - t39 * t56;
	t29 = sin(pkin(7));
	t32 = cos(pkin(7));
	t30 = sin(pkin(6));
	t69 = t30 * t40;
	t15 = -t21 * t29 + t32 * t69;
	t28 = sin(pkin(8));
	t31 = cos(pkin(8));
	t22 = t35 * t56 + t36 * t39;
	t34 = sin(qJ(3));
	t38 = cos(qJ(3));
	t59 = t29 * t69;
	t5 = (t21 * t32 + t59) * t38 + t22 * t34;
	t52 = t15 * t28 + t31 * t5;
	t66 = t32 * t34;
	t6 = t21 * t66 - t22 * t38 + t34 * t59;
	t82 = t52 * t33 + t6 * t37;
	t81 = -t6 * t33 + t52 * t37;
	t54 = t33 * r_i_i_C(1) + t37 * r_i_i_C(2);
	t75 = pkin(12) + r_i_i_C(3);
	t58 = t75 * t28;
	t78 = -t54 * t31 + t58;
	t74 = t29 * pkin(11);
	t72 = t28 * t29;
	t71 = t29 * t31;
	t70 = t30 * t36;
	t68 = t31 * t33;
	t67 = t31 * t37;
	t65 = t32 * t38;
	t64 = t34 * t35;
	t63 = t34 * t39;
	t62 = t35 * t38;
	t61 = t38 * t39;
	t57 = t36 * t60;
	t55 = t60 * t29;
	t23 = -t40 * t35 - t39 * t57;
	t17 = -t23 * t29 + t32 * t70;
	t24 = -t35 * t57 + t40 * t39;
	t41 = t23 * t32 + t29 * t70;
	t7 = -t24 * t34 + t41 * t38;
	t49 = t17 * t28 + t31 * t7;
	t13 = t38 * t55 + (t32 * t61 - t64) * t30;
	t48 = t13 * t31 + (-t30 * t39 * t29 + t60 * t32) * t28;
	t47 = t37 * r_i_i_C(1) - t33 * r_i_i_C(2) + pkin(3);
	t9 = t21 * t34 - t22 * t65;
	t45 = t22 * t72 + t31 * t9;
	t11 = -t23 * t34 - t24 * t65;
	t43 = t11 * t31 + t24 * t72;
	t14 = t34 * t55 + (t32 * t63 + t62) * t30;
	t12 = t23 * t38 - t24 * t66;
	t10 = -t21 * t38 - t22 * t66;
	t8 = t24 * t38 + t41 * t34;
	t2 = t49 * t33 + t8 * t37;
	t1 = -t8 * t33 + t49 * t37;
	t3 = [t82 * r_i_i_C(1) + t81 * r_i_i_C(2) + t6 * pkin(3) - t22 * pkin(2) - t36 * pkin(1) + pkin(10) * t69 + t15 * pkin(11) + t75 * (t15 * t31 - t5 * t28), (t12 * t37 + t43 * t33) * r_i_i_C(1) + (-t12 * t33 + t43 * t37) * r_i_i_C(2) + t12 * pkin(3) + t23 * pkin(2) + t24 * t74 + t75 * (-t11 * t28 + t24 * t71), (t7 * t37 - t8 * t68) * r_i_i_C(1) + (-t7 * t33 - t8 * t67) * r_i_i_C(2) + t7 * pkin(3) + t8 * t58, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t40 * pkin(1) + t24 * pkin(2) + t8 * pkin(3) + pkin(10) * t70 + t17 * pkin(11) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t75 * (t17 * t31 - t7 * t28), (t10 * t37 + t45 * t33) * r_i_i_C(1) + (-t10 * t33 + t45 * t37) * r_i_i_C(2) + t10 * pkin(3) - t21 * pkin(2) + t22 * t74 + t75 * (t22 * t71 - t9 * t28), (-t37 * t5 + t6 * t68) * r_i_i_C(1) + (t33 * t5 + t6 * t67) * r_i_i_C(2) - t5 * pkin(3) - t6 * t58, -t81 * r_i_i_C(1) + t82 * r_i_i_C(2), 0, 0; 0, (t47 * (-t32 * t64 + t61) - t78 * (-t32 * t62 - t63) + t39 * pkin(2) + (t54 * t28 + t75 * t31 + pkin(11)) * t35 * t29) * t30, t47 * t13 + t78 * t14, (-t14 * t33 + t48 * t37) * r_i_i_C(1) + (-t14 * t37 - t48 * t33) * r_i_i_C(2), 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:54
	% DurationCPUTime: 1.26s
	% Computational Cost: add. (840->137), mult. (2410->242), div. (0->0), fcn. (3166->16), ass. (0->89)
	t64 = sin(qJ(2));
	t65 = sin(qJ(1));
	t69 = cos(qJ(2));
	t70 = cos(qJ(1));
	t87 = cos(pkin(6));
	t83 = t70 * t87;
	t49 = t65 * t64 - t69 * t83;
	t50 = t64 * t83 + t65 * t69;
	t60 = cos(pkin(7));
	t63 = sin(qJ(3));
	t68 = cos(qJ(3));
	t57 = sin(pkin(7));
	t58 = sin(pkin(6));
	t96 = t58 * t70;
	t85 = t57 * t96;
	t32 = (t49 * t60 + t85) * t68 + t50 * t63;
	t44 = -t49 * t57 + t60 * t96;
	t56 = sin(pkin(8));
	t59 = cos(pkin(8));
	t19 = t32 * t56 - t44 * t59;
	t93 = t60 * t63;
	t33 = t49 * t93 - t50 * t68 + t63 * t85;
	t62 = sin(qJ(4));
	t67 = cos(qJ(4));
	t80 = t32 * t59 + t44 * t56;
	t6 = t33 * t67 + t80 * t62;
	t61 = sin(qJ(5));
	t66 = cos(qJ(5));
	t116 = t19 * t66 + t6 * t61;
	t115 = -t19 * t61 + t6 * t66;
	t112 = t33 * t62 - t80 * t67;
	t107 = r_i_i_C(3) + pkin(13);
	t106 = t56 * pkin(12);
	t105 = t57 * pkin(11);
	t89 = t64 * t68;
	t90 = t63 * t69;
	t46 = (-t60 * t89 - t90) * t58;
	t104 = t46 * t56;
	t102 = t56 * t57;
	t101 = t56 * t61;
	t100 = t56 * t66;
	t99 = t57 * t59;
	t98 = t57 * t64;
	t97 = t58 * t65;
	t95 = t59 * t62;
	t94 = t59 * t67;
	t92 = t60 * t68;
	t91 = t63 * t64;
	t88 = t68 * t69;
	t86 = t58 * t98;
	t84 = t65 * t87;
	t82 = t87 * t57;
	t42 = t68 * t82 + (t60 * t88 - t91) * t58;
	t48 = -t58 * t69 * t57 + t87 * t60;
	t79 = t42 * t59 + t48 * t56;
	t78 = t66 * r_i_i_C(1) - t61 * r_i_i_C(2) + pkin(4);
	t36 = t49 * t63 - t50 * t92;
	t26 = -t36 * t56 + t50 * t99;
	t77 = t50 * t102 + t36 * t59;
	t51 = -t70 * t64 - t69 * t84;
	t52 = -t64 * t84 + t70 * t69;
	t38 = -t51 * t63 - t52 * t92;
	t27 = -t38 * t56 + t52 * t99;
	t76 = t52 * t102 + t38 * t59;
	t74 = -t51 * t57 + t60 * t97;
	t73 = t51 * t60 + t57 * t97;
	t72 = t46 * t59 + t56 * t86;
	t71 = t74 * t56;
	t34 = -t52 * t63 + t73 * t68;
	t21 = -t34 * t56 + t74 * t59;
	t47 = (-t60 * t91 + t88) * t58;
	t43 = t63 * t82 + (t60 * t90 + t89) * t58;
	t40 = t59 * t86 - t104;
	t39 = t51 * t68 - t52 * t93;
	t37 = -t49 * t68 - t50 * t93;
	t35 = t52 * t68 + t73 * t63;
	t29 = -t42 * t56 + t48 * t59;
	t25 = t47 * t67 + t72 * t62;
	t23 = t42 * t67 - t43 * t95;
	t18 = t43 * t67 + t79 * t62;
	t16 = t34 * t67 - t35 * t95;
	t14 = -t32 * t67 + t33 * t95;
	t12 = t39 * t67 + t76 * t62;
	t10 = t37 * t67 + t77 * t62;
	t8 = t35 * t67 + (t34 * t59 + t71) * t62;
	t7 = -t34 * t94 + t35 * t62 - t67 * t71;
	t2 = t21 * t61 + t8 * t66;
	t1 = t21 * t66 - t8 * t61;
	t3 = [-t65 * pkin(1) - t50 * pkin(2) + t33 * pkin(3) + t6 * pkin(4) + pkin(10) * t96 + t44 * pkin(11) - t19 * pkin(12) + t115 * r_i_i_C(1) - t116 * r_i_i_C(2) + t107 * t112, (t12 * t66 + t27 * t61) * r_i_i_C(1) + (-t12 * t61 + t27 * t66) * r_i_i_C(2) + t12 * pkin(4) + t39 * pkin(3) + t51 * pkin(2) + t52 * t105 + t107 * (t39 * t62 - t76 * t67) + t27 * pkin(12), (t35 * t101 + t16 * t66) * r_i_i_C(1) + (t35 * t100 - t16 * t61) * r_i_i_C(2) + t16 * pkin(4) + t34 * pkin(3) + t35 * t106 + t107 * (t34 * t62 + t35 * t94), t107 * t8 - t78 * t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t70 * pkin(1) + t52 * pkin(2) + t35 * pkin(3) + t8 * pkin(4) + pkin(10) * t97 + t74 * pkin(11) + t21 * pkin(12) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t107 * t7, (t10 * t66 + t26 * t61) * r_i_i_C(1) + (-t10 * t61 + t26 * t66) * r_i_i_C(2) + t10 * pkin(4) + t37 * pkin(3) - t49 * pkin(2) + t50 * t105 + t107 * (t37 * t62 - t77 * t67) + t26 * pkin(12), (-t101 * t33 + t14 * t66) * r_i_i_C(1) + (-t100 * t33 - t14 * t61) * r_i_i_C(2) + t14 * pkin(4) - t32 * pkin(3) - t33 * t106 + t107 * (-t32 * t62 - t33 * t94), -t107 * t6 + t112 * t78, t116 * r_i_i_C(1) + t115 * r_i_i_C(2), 0; 0, (t25 * t66 + t40 * t61) * r_i_i_C(1) + (-t25 * t61 + t40 * t66) * r_i_i_C(2) + t25 * pkin(4) + t47 * pkin(3) - pkin(12) * t104 + t107 * (t47 * t62 - t72 * t67) + (t69 * pkin(2) + (pkin(12) * t59 + pkin(11)) * t98) * t58, (t43 * t101 + t23 * t66) * r_i_i_C(1) + (t43 * t100 - t23 * t61) * r_i_i_C(2) + t23 * pkin(4) + t42 * pkin(3) + t43 * t106 + t107 * (t42 * t62 + t43 * t94), t107 * t18 + t78 * (-t43 * t62 + t79 * t67), (-t18 * t61 + t29 * t66) * r_i_i_C(1) + (-t18 * t66 - t29 * t61) * r_i_i_C(2), 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:55
	% DurationCPUTime: 1.96s
	% Computational Cost: add. (1756->187), mult. (5055->316), div. (0->0), fcn. (6693->18), ass. (0->115)
	t76 = sin(pkin(6));
	t88 = cos(qJ(1));
	t128 = t76 * t88;
	t120 = cos(pkin(6));
	t113 = t88 * t120;
	t138 = sin(qJ(1));
	t83 = sin(qJ(2));
	t87 = cos(qJ(2));
	t68 = -t87 * t113 + t138 * t83;
	t75 = sin(pkin(7));
	t78 = cos(pkin(7));
	t64 = t78 * t128 - t68 * t75;
	t74 = sin(pkin(8));
	t137 = t64 * t74;
	t139 = cos(qJ(4));
	t118 = t75 * t128;
	t69 = t83 * t113 + t138 * t87;
	t82 = sin(qJ(3));
	t86 = cos(qJ(3));
	t143 = (t68 * t78 + t118) * t86 + t69 * t82;
	t126 = t78 * t82;
	t53 = t82 * t118 + t68 * t126 - t69 * t86;
	t77 = cos(pkin(8));
	t81 = sin(qJ(4));
	t22 = t53 * t139 + (t143 * t77 + t137) * t81;
	t40 = t143 * t74 - t64 * t77;
	t80 = sin(qJ(5));
	t85 = cos(qJ(5));
	t6 = t22 * t85 - t40 * t80;
	t79 = sin(qJ(6));
	t151 = t6 * t79;
	t84 = cos(qJ(6));
	t150 = t6 * t84;
	t149 = t22 * t80 + t40 * t85;
	t146 = t53 * t81;
	t142 = r_i_i_C(3) + pkin(14);
	t110 = t120 * t138;
	t70 = -t83 * t110 + t87 * t88;
	t104 = t87 * t110 + t88 * t83;
	t116 = t76 * t138;
	t96 = -t104 * t78 + t75 * t116;
	t92 = t70 * t82 - t96 * t86;
	t97 = t104 * t75 + t78 * t116;
	t89 = t92 * t74 + t97 * t77;
	t141 = pkin(11) * t75;
	t140 = pkin(12) * t74;
	t122 = t83 * t86;
	t123 = t82 * t87;
	t66 = (-t78 * t122 - t123) * t76;
	t135 = t66 * t74;
	t133 = t74 * t75;
	t132 = t74 * t80;
	t131 = t74 * t85;
	t130 = t75 * t77;
	t129 = t75 * t83;
	t127 = t77 * t81;
	t125 = t78 * t86;
	t124 = t82 * t83;
	t121 = t86 * t87;
	t119 = t76 * t129;
	t117 = t74 * t139;
	t115 = t77 * t139;
	t114 = t75 * t120;
	t112 = t75 * t117;
	t109 = r_i_i_C(1) * t84 - r_i_i_C(2) * t79 + pkin(5);
	t108 = r_i_i_C(1) * t79 + r_i_i_C(2) * t84 + pkin(13);
	t55 = -t69 * t125 + t68 * t82;
	t46 = t69 * t130 - t55 * t74;
	t57 = t104 * t82 - t70 * t125;
	t47 = t70 * t130 - t57 * t74;
	t105 = -t76 * t87 * t75 + t120 * t78;
	t103 = t105 * t74;
	t100 = -t109 * t85 - t142 * t80 - pkin(4);
	t99 = t86 * t114 + (t78 * t121 - t124) * t76;
	t98 = t143 * t139;
	t95 = t99 * t139;
	t94 = t97 * t74;
	t90 = t92 * t139;
	t67 = (-t78 * t124 + t121) * t76;
	t63 = t82 * t114 + (t78 * t123 + t122) * t76;
	t59 = t77 * t119 - t135;
	t58 = -t104 * t86 - t70 * t126;
	t56 = -t69 * t126 - t68 * t86;
	t54 = t70 * t86 + t96 * t82;
	t50 = t105 * t77 - t99 * t74;
	t45 = t67 * t139 + (t74 * t119 + t66 * t77) * t81;
	t44 = -t76 * t83 * t112 - t66 * t115 + t67 * t81;
	t43 = -t63 * t127 + t95;
	t42 = t63 * t115 + t99 * t81;
	t38 = t63 * t139 + (t99 * t77 + t103) * t81;
	t37 = -t139 * t103 + t63 * t81 - t77 * t95;
	t36 = t63 * t132 + t43 * t85;
	t34 = t45 * t85 + t59 * t80;
	t32 = -t54 * t127 - t90;
	t31 = t54 * t115 - t92 * t81;
	t30 = t127 * t53 - t98;
	t29 = -t115 * t53 - t143 * t81;
	t28 = t58 * t139 + (t70 * t133 + t57 * t77) * t81;
	t27 = -t70 * t112 - t57 * t115 + t58 * t81;
	t26 = t56 * t139 + (t69 * t133 + t55 * t77) * t81;
	t25 = -t69 * t112 - t55 * t115 + t56 * t81;
	t24 = t54 * t139 + (-t92 * t77 + t94) * t81;
	t23 = -t139 * t94 + t54 * t81 + t77 * t90;
	t21 = -t115 * t143 - t64 * t117 + t146;
	t19 = t137 * t139 + t77 * t98 - t146;
	t18 = t38 * t85 + t50 * t80;
	t16 = t54 * t132 + t32 * t85;
	t14 = -t132 * t53 + t30 * t85;
	t12 = t28 * t85 + t47 * t80;
	t10 = t26 * t85 + t46 * t80;
	t8 = t24 * t85 + t89 * t80;
	t7 = t24 * t80 - t89 * t85;
	t2 = t23 * t79 + t8 * t84;
	t1 = t23 * t84 - t79 * t8;
	t3 = [(t21 * t79 + t150) * r_i_i_C(1) + (t21 * t84 - t151) * r_i_i_C(2) + t6 * pkin(5) + t22 * pkin(4) + t21 * pkin(13) + t53 * pkin(3) - t69 * pkin(2) - t138 * pkin(1) + pkin(10) * t128 + t142 * t149 - t40 * pkin(12) + t64 * pkin(11), (t12 * t84 + t27 * t79) * r_i_i_C(1) + (-t12 * t79 + t27 * t84) * r_i_i_C(2) + t12 * pkin(5) + t28 * pkin(4) + t27 * pkin(13) + t58 * pkin(3) - t104 * pkin(2) + t70 * t141 + t142 * (t28 * t80 - t47 * t85) + t47 * pkin(12), (t16 * t84 + t31 * t79) * r_i_i_C(1) + (-t16 * t79 + t31 * t84) * r_i_i_C(2) + t16 * pkin(5) + t32 * pkin(4) + t31 * pkin(13) - t92 * pkin(3) + t54 * t140 + t142 * (-t54 * t131 + t32 * t80), t100 * t23 + t108 * t24, -t109 * t7 + t142 * t8, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t88 * pkin(1) + t70 * pkin(2) + t54 * pkin(3) + t24 * pkin(4) + t8 * pkin(5) + pkin(10) * t116 + t97 * pkin(11) + t89 * pkin(12) + t23 * pkin(13) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t142 * t7, (t10 * t84 + t25 * t79) * r_i_i_C(1) + (-t10 * t79 + t25 * t84) * r_i_i_C(2) + t10 * pkin(5) + t26 * pkin(4) + t25 * pkin(13) + t56 * pkin(3) - t68 * pkin(2) + t69 * t141 + t142 * (t26 * t80 - t46 * t85) + t46 * pkin(12), (t14 * t84 + t29 * t79) * r_i_i_C(1) + (-t14 * t79 + t29 * t84) * r_i_i_C(2) + t14 * pkin(5) + t30 * pkin(4) + t29 * pkin(13) - t143 * pkin(3) - t53 * t140 + t142 * (t131 * t53 + t30 * t80), t100 * t19 - t108 * t22, t109 * t149 - t142 * t6, (t19 * t84 + t151) * r_i_i_C(1) + (-t19 * t79 + t150) * r_i_i_C(2); 0, (t34 * t84 + t44 * t79) * r_i_i_C(1) + (-t34 * t79 + t44 * t84) * r_i_i_C(2) + t34 * pkin(5) + t45 * pkin(4) + t44 * pkin(13) + t67 * pkin(3) - pkin(12) * t135 + t142 * (t45 * t80 - t59 * t85) + (t87 * pkin(2) + (pkin(12) * t77 + pkin(11)) * t129) * t76, (t36 * t84 + t42 * t79) * r_i_i_C(1) + (-t36 * t79 + t42 * t84) * r_i_i_C(2) + t36 * pkin(5) + t43 * pkin(4) + t42 * pkin(13) + t99 * pkin(3) + t63 * t140 + t142 * (-t63 * t131 + t43 * t80), t100 * t37 + t108 * t38, t142 * t18 + t109 * (-t38 * t80 + t50 * t85), (-t18 * t79 + t37 * t84) * r_i_i_C(1) + (-t18 * t84 - t37 * t79) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end