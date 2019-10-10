% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRR8
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:29
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR8_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR8_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR8_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobia_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
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
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
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
	t11 = (pkin(9) + r_i_i_C(3)) * t5;
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
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (93->46), mult. (252->84), div. (0->0), fcn. (316->10), ass. (0->33)
	t39 = pkin(10) + r_i_i_C(3);
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
	t7 = [t19 * r_i_i_C(1) + t38 * r_i_i_C(2) - t4 * pkin(2) - t15 * pkin(1) + pkin(9) * t35 + t39 * (-t3 * t10 + t12 * t35), (t5 * t16 - t6 * t34) * r_i_i_C(1) + (-t5 * t13 - t6 * t33) * r_i_i_C(2) + t5 * pkin(2) + t6 * t26, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; t18 * pkin(1) + t6 * pkin(2) + pkin(9) * t36 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t39 * (-t5 * t10 + t12 * t36), (-t3 * t16 - t4 * t34) * r_i_i_C(1) + (t3 * t13 - t4 * t33) * r_i_i_C(2) - t3 * pkin(2) + t4 * t26, -t38 * r_i_i_C(1) + t19 * r_i_i_C(2), 0, 0, 0; 0, ((-t12 * t32 + t29) * r_i_i_C(1) + (-t12 * t30 - t31) * r_i_i_C(2) + t17 * pkin(2) + t14 * t26) * t11, ((t12 * t29 - t32) * r_i_i_C(1) + (-t12 * t31 - t30) * r_i_i_C(2)) * t11 + (t16 * r_i_i_C(1) - t13 * r_i_i_C(2)) * t10 * t28, 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:31
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (242->73), mult. (663->129), div. (0->0), fcn. (854->12), ass. (0->49)
	t35 = sin(qJ(2));
	t36 = sin(qJ(1));
	t39 = cos(qJ(2));
	t40 = cos(qJ(1));
	t48 = cos(pkin(6));
	t44 = t40 * t48;
	t22 = t36 * t35 - t39 * t44;
	t30 = sin(pkin(7));
	t32 = cos(pkin(7));
	t31 = sin(pkin(6));
	t55 = t31 * t40;
	t15 = -t22 * t30 + t32 * t55;
	t33 = sin(qJ(4));
	t37 = cos(qJ(4));
	t23 = t35 * t44 + t36 * t39;
	t34 = sin(qJ(3));
	t38 = cos(qJ(3));
	t46 = t30 * t55;
	t54 = t32 * t34;
	t6 = t22 * t54 - t23 * t38 + t34 * t46;
	t62 = t15 * t37 - t6 * t33;
	t61 = t15 * t33 + t6 * t37;
	t60 = r_i_i_C(3) + pkin(11);
	t59 = t30 * pkin(10);
	t58 = t30 * t33;
	t57 = t30 * t37;
	t56 = t31 * t36;
	t53 = t32 * t38;
	t52 = t34 * t35;
	t51 = t34 * t39;
	t50 = t35 * t38;
	t49 = t38 * t39;
	t47 = t30 * t56;
	t45 = t36 * t48;
	t43 = t48 * t30;
	t42 = t37 * r_i_i_C(1) - t33 * r_i_i_C(2) + pkin(3);
	t24 = -t40 * t35 - t39 * t45;
	t17 = -t24 * t30 + t32 * t56;
	t41 = -t23 * t34 + (-t22 * t32 - t46) * t38;
	t25 = -t35 * t45 + t40 * t39;
	t21 = -t31 * t39 * t30 + t48 * t32;
	t14 = t34 * t43 + (t32 * t51 + t50) * t31;
	t12 = t24 * t38 - t25 * t54;
	t10 = -t22 * t38 - t23 * t54;
	t8 = t25 * t38 + (t24 * t32 + t47) * t34;
	t7 = -t24 * t53 + t25 * t34 - t38 * t47;
	t2 = t17 * t33 + t8 * t37;
	t1 = t17 * t37 - t8 * t33;
	t3 = [-t36 * pkin(1) - t23 * pkin(2) + t6 * pkin(3) + pkin(9) * t55 + t15 * pkin(10) + t61 * r_i_i_C(1) + t62 * r_i_i_C(2) + t60 * t41, (t12 * t37 + t25 * t58) * r_i_i_C(1) + (-t12 * t33 + t25 * t57) * r_i_i_C(2) + t12 * pkin(3) + t24 * pkin(2) + t25 * t59 + t60 * (t24 * t34 + t25 * t53), -t42 * t7 + t60 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t40 * pkin(1) + t25 * pkin(2) + t8 * pkin(3) + pkin(9) * t56 + t17 * pkin(10) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t60 * t7, (t10 * t37 + t23 * t58) * r_i_i_C(1) + (-t10 * t33 + t23 * t57) * r_i_i_C(2) + t10 * pkin(3) - t22 * pkin(2) + t23 * t59 + t60 * (-t22 * t34 + t23 * t53), t42 * t41 - t6 * t60, -t62 * r_i_i_C(1) + t61 * r_i_i_C(2), 0, 0; 0, (t42 * (-t32 * t52 + t49) + t39 * pkin(2) + (t33 * r_i_i_C(1) + t37 * r_i_i_C(2) + pkin(10)) * t35 * t30 + t60 * (t32 * t50 + t51)) * t31, t60 * t14 + t42 * (t38 * t43 + (t32 * t49 - t52) * t31), (-t14 * t33 + t21 * t37) * r_i_i_C(1) + (-t14 * t37 - t21 * t33) * r_i_i_C(2), 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (379->72), mult. (871->123), div. (0->0), fcn. (1119->14), ass. (0->52)
	t43 = sin(qJ(4));
	t75 = pkin(4) * t43 + pkin(10);
	t47 = cos(qJ(4));
	t36 = t47 * pkin(4) + pkin(3);
	t39 = qJ(4) + qJ(5);
	t37 = sin(t39);
	t38 = cos(t39);
	t54 = t38 * r_i_i_C(1) - t37 * r_i_i_C(2) + t36;
	t74 = r_i_i_C(1) * t37 + r_i_i_C(2) * t38 + t75;
	t45 = sin(qJ(2));
	t46 = sin(qJ(1));
	t49 = cos(qJ(2));
	t50 = cos(qJ(1));
	t60 = cos(pkin(6));
	t56 = t50 * t60;
	t28 = t46 * t45 - t49 * t56;
	t29 = t45 * t56 + t46 * t49;
	t44 = sin(qJ(3));
	t48 = cos(qJ(3));
	t40 = sin(pkin(7));
	t41 = sin(pkin(6));
	t67 = t41 * t50;
	t58 = t40 * t67;
	t42 = cos(pkin(7));
	t66 = t42 * t44;
	t12 = t28 * t66 - t29 * t48 + t44 * t58;
	t21 = -t28 * t40 + t42 * t67;
	t73 = (t12 * t37 - t21 * t38) * r_i_i_C(1) + (t12 * t38 + t21 * t37) * r_i_i_C(2);
	t57 = t46 * t60;
	t30 = -t50 * t45 - t49 * t57;
	t31 = -t45 * t57 + t50 * t49;
	t68 = t41 * t46;
	t59 = t40 * t68;
	t14 = t31 * t48 + (t30 * t42 + t59) * t44;
	t23 = -t30 * t40 + t42 * t68;
	t5 = -t14 * t37 + t23 * t38;
	t6 = t14 * t38 + t23 * t37;
	t72 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t55 = t60 * t40;
	t62 = t45 * t48;
	t63 = t44 * t49;
	t20 = t44 * t55 + (t42 * t63 + t62) * t41;
	t27 = -t41 * t49 * t40 + t60 * t42;
	t71 = (-t20 * t37 + t27 * t38) * r_i_i_C(1) + (-t20 * t38 - t27 * t37) * r_i_i_C(2);
	t69 = r_i_i_C(3) + pkin(12) + pkin(11);
	t65 = t42 * t48;
	t64 = t44 * t45;
	t61 = t48 * t49;
	t53 = t74 * t40;
	t52 = -t29 * t44 + (-t28 * t42 - t58) * t48;
	t13 = -t30 * t65 + t31 * t44 - t48 * t59;
	t1 = [-t46 * pkin(1) - t29 * pkin(2) + pkin(9) * t67 + t54 * t12 + t74 * t21 + t69 * t52, t30 * pkin(2) + t54 * (t30 * t48 - t31 * t66) + t31 * t53 + t69 * (t30 * t44 + t31 * t65), -t54 * t13 + t69 * t14, (-t14 * t43 + t23 * t47) * pkin(4) + t72, t72, 0; t50 * pkin(1) + t31 * pkin(2) + pkin(9) * t68 + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t69 * t13 + t14 * t36 + t75 * t23, -t28 * pkin(2) + t54 * (-t28 * t48 - t29 * t66) + t29 * t53 + t69 * (-t28 * t44 + t29 * t65), -t12 * t69 + t54 * t52, (t12 * t43 - t21 * t47) * pkin(4) + t73, t73, 0; 0, (t69 * (t42 * t62 + t63) + t54 * (-t42 * t64 + t61) + t49 * pkin(2) + t45 * t53) * t41, t69 * t20 + t54 * (t48 * t55 + (t42 * t61 - t64) * t41), (-t20 * t43 + t27 * t47) * pkin(4) + t71, t71, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:31
	% DurationCPUTime: 0.80s
	% Computational Cost: add. (822->114), mult. (1862->190), div. (0->0), fcn. (2429->16), ass. (0->78)
	t64 = sin(qJ(2));
	t67 = cos(qJ(2));
	t68 = cos(qJ(1));
	t89 = cos(pkin(6));
	t82 = t68 * t89;
	t97 = sin(qJ(1));
	t46 = t97 * t64 - t67 * t82;
	t47 = t64 * t82 + t97 * t67;
	t63 = sin(qJ(3));
	t58 = sin(pkin(7));
	t59 = sin(pkin(6));
	t93 = t59 * t68;
	t87 = t58 * t93;
	t60 = cos(pkin(7));
	t92 = t60 * t63;
	t98 = cos(qJ(3));
	t25 = -t46 * t92 + t47 * t98 - t63 * t87;
	t39 = -t46 * t58 + t60 * t93;
	t57 = qJ(4) + qJ(5);
	t55 = sin(t57);
	t56 = cos(t57);
	t10 = t25 * t56 - t39 * t55;
	t85 = t60 * t98;
	t24 = t46 * t85 + t47 * t63 + t98 * t87;
	t61 = sin(qJ(6));
	t65 = cos(qJ(6));
	t109 = t10 * t61 - t24 * t65;
	t108 = -t10 * t65 - t24 * t61;
	t102 = r_i_i_C(3) + pkin(13);
	t105 = t65 * r_i_i_C(1) - t61 * r_i_i_C(2) + pkin(5);
	t62 = sin(qJ(4));
	t104 = pkin(4) * t62 + pkin(10);
	t66 = cos(qJ(4));
	t54 = t66 * pkin(4) + pkin(3);
	t103 = t102 * t55 + t54;
	t96 = t55 * t58;
	t95 = t56 * t58;
	t94 = t58 * t59;
	t91 = t63 * t64;
	t90 = t63 * t67;
	t88 = t64 * t94;
	t86 = t59 * t97;
	t84 = t98 * t64;
	t83 = t98 * t67;
	t81 = t89 * t58;
	t80 = t58 * t86;
	t79 = t104 * t58;
	t78 = t89 * t97;
	t69 = -pkin(12) - pkin(11);
	t77 = t61 * r_i_i_C(1) + t65 * r_i_i_C(2) - t69;
	t76 = t102 * t10 + t105 * (-t25 * t55 - t39 * t56);
	t48 = -t64 * t78 + t68 * t67;
	t73 = t68 * t64 + t67 * t78;
	t29 = t48 * t98 + (-t73 * t60 + t80) * t63;
	t70 = t73 * t58 + t60 * t86;
	t13 = t29 * t55 - t70 * t56;
	t14 = t29 * t56 + t70 * t55;
	t75 = t102 * t14 - t105 * t13;
	t38 = t63 * t81 + (t60 * t90 + t84) * t59;
	t45 = t89 * t60 - t67 * t94;
	t23 = t38 * t56 + t45 * t55;
	t74 = t102 * t23 + t105 * (-t38 * t55 + t45 * t56);
	t72 = t73 * t98;
	t71 = -t105 * t56 - t103;
	t44 = (-t60 * t91 + t83) * t59;
	t43 = (t60 * t84 + t90) * t59;
	t37 = -t98 * t81 + (-t60 * t83 + t91) * t59;
	t35 = t44 * t56 + t55 * t88;
	t33 = -t48 * t92 - t72;
	t32 = t48 * t85 - t73 * t63;
	t31 = -t46 * t98 - t47 * t92;
	t30 = -t46 * t63 + t47 * t85;
	t28 = t48 * t63 + t60 * t72 - t98 * t80;
	t18 = t33 * t56 + t48 * t96;
	t16 = t31 * t56 + t47 * t96;
	t2 = t14 * t65 + t28 * t61;
	t1 = -t14 * t61 + t28 * t65;
	t3 = [t108 * r_i_i_C(1) + t109 * r_i_i_C(2) - t10 * pkin(5) + t24 * t69 - t47 * pkin(2) - t97 * pkin(1) + pkin(9) * t93 - t103 * t25 + (-t102 * t56 + t104) * t39, (t18 * t65 + t32 * t61) * r_i_i_C(1) + (-t18 * t61 + t32 * t65) * r_i_i_C(2) + t18 * pkin(5) + t33 * t54 - t32 * t69 - t73 * pkin(2) + t48 * t79 + t102 * (t33 * t55 - t48 * t95), t71 * t28 + t77 * t29, (-t29 * t62 + t70 * t66) * pkin(4) + t75, t75, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t68 * pkin(1) + t48 * pkin(2) + t14 * pkin(5) + pkin(9) * t86 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t102 * t13 + t104 * t70 - t28 * t69 + t29 * t54, (t16 * t65 + t30 * t61) * r_i_i_C(1) + (-t16 * t61 + t30 * t65) * r_i_i_C(2) + t16 * pkin(5) + t31 * t54 - t30 * t69 - t46 * pkin(2) + t47 * t79 + t102 * (t31 * t55 - t47 * t95), t71 * t24 + t77 * t25, (-t25 * t62 - t39 * t66) * pkin(4) + t76, t76, -t109 * r_i_i_C(1) + t108 * r_i_i_C(2); 0, (t35 * t65 + t43 * t61) * r_i_i_C(1) + (-t35 * t61 + t43 * t65) * r_i_i_C(2) + t35 * pkin(5) + t44 * t54 - t43 * t69 + t102 * (t44 * t55 - t56 * t88) + (t67 * pkin(2) + t64 * t79) * t59, t71 * t37 + t77 * t38, (-t38 * t62 + t45 * t66) * pkin(4) + t74, t74, (-t23 * t61 + t37 * t65) * r_i_i_C(1) + (-t23 * t65 - t37 * t61) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end