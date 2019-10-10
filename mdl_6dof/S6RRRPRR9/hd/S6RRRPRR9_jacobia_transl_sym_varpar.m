% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:08
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR9_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR9_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR9_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobia_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
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
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
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
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.29s
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
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (172->58), mult. (443->100), div. (0->0), fcn. (561->12), ass. (0->39)
	t24 = sin(qJ(3));
	t43 = t24 * pkin(3);
	t27 = cos(qJ(3));
	t42 = t27 * pkin(3);
	t21 = sin(pkin(6));
	t26 = sin(qJ(1));
	t41 = t21 * t26;
	t28 = cos(qJ(2));
	t40 = t21 * t28;
	t29 = cos(qJ(1));
	t39 = t21 * t29;
	t38 = pkin(10) + qJ(4);
	t37 = cos(pkin(6));
	t36 = t26 * t37;
	t35 = t29 * t37;
	t20 = sin(pkin(7));
	t23 = cos(pkin(7));
	t19 = sin(pkin(13));
	t22 = cos(pkin(13));
	t32 = t27 * t19 + t24 * t22;
	t4 = t32 * t20;
	t34 = r_i_i_C(1) * t4 + t20 * t43 + t38 * t23 + pkin(9);
	t25 = sin(qJ(2));
	t11 = -t29 * t25 - t28 * t36;
	t12 = -t25 * t36 + t29 * t28;
	t13 = t24 * t19 - t27 * t22;
	t6 = t32 * t23;
	t33 = t11 * t6 - t12 * t13;
	t18 = pkin(2) + t42;
	t31 = t13 * r_i_i_C(1) + r_i_i_C(2) * t32 - t18;
	t5 = t13 * t23;
	t8 = -t38 * t20 + t23 * t43;
	t30 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2) - t20 * r_i_i_C(3) + t8;
	t10 = t25 * t35 + t26 * t28;
	t9 = t26 * t25 - t28 * t35;
	t3 = t13 * t20;
	t2 = -t11 * t20 + t23 * t41;
	t1 = -t11 * t5 - t12 * t32 - t3 * t41;
	t7 = [-t26 * pkin(1) + t31 * t10 + t30 * t9 + (-r_i_i_C(2) * t3 + r_i_i_C(3) * t23 + t34) * t39, -t31 * t11 - t30 * t12, t1 * r_i_i_C(1) + (-t4 * t41 - t33) * r_i_i_C(2) + (-t12 * t24 + (t11 * t23 + t20 * t41) * t27) * pkin(3), t2, 0, 0; t29 * pkin(1) + t33 * r_i_i_C(1) + t1 * r_i_i_C(2) + t2 * r_i_i_C(3) + t11 * t8 + t12 * t18 + t34 * t41, -t30 * t10 + t31 * t9, (-t10 * t32 + t3 * t39 + t9 * t5) * r_i_i_C(1) + (t10 * t13 + t4 * t39 + t9 * t6) * r_i_i_C(2) + (-t10 * t24 + (-t20 * t39 - t23 * t9) * t27) * pkin(3), t9 * t20 - t23 * t39, 0, 0; 0, (-t30 * t25 - t31 * t28) * t21, (-t37 * t3 + (-t25 * t32 - t28 * t5) * t21) * r_i_i_C(1) + (-t37 * t4 + (t13 * t25 - t28 * t6) * t21) * r_i_i_C(2) - t21 * t25 * t43 + (t37 * t20 + t23 * t40) * t42, -t20 * t40 + t37 * t23, 0, 0;];
	Ja_transl = t7;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.48s
	% Computational Cost: add. (402->90), mult. (1067->155), div. (0->0), fcn. (1393->14), ass. (0->52)
	t46 = sin(qJ(2));
	t47 = sin(qJ(1));
	t50 = cos(qJ(2));
	t51 = cos(qJ(1));
	t61 = cos(pkin(6));
	t57 = t51 * t61;
	t30 = t47 * t46 - t50 * t57;
	t41 = sin(pkin(7));
	t43 = cos(pkin(7));
	t42 = sin(pkin(6));
	t63 = t42 * t51;
	t20 = -t30 * t41 + t43 * t63;
	t44 = sin(qJ(5));
	t48 = cos(qJ(5));
	t40 = sin(pkin(13));
	t45 = sin(qJ(3));
	t49 = cos(qJ(3));
	t60 = cos(pkin(13));
	t35 = -t49 * t40 - t45 * t60;
	t24 = t35 * t41;
	t26 = t35 * t43;
	t31 = t46 * t57 + t47 * t50;
	t53 = -t45 * t40 + t49 * t60;
	t7 = -t24 * t63 - t26 * t30 - t31 * t53;
	t70 = t20 * t48 - t44 * t7;
	t69 = t20 * t44 + t48 * t7;
	t13 = (-t26 * t50 + t46 * t53) * t42 - t61 * t24;
	t68 = -r_i_i_C(3) - pkin(11);
	t67 = pkin(3) * t45;
	t66 = t41 * t44;
	t65 = t41 * t48;
	t64 = t42 * t47;
	t62 = pkin(10) + qJ(4);
	t59 = t42 * (t41 * t67 + t43 * t62 + pkin(9));
	t58 = t47 * t61;
	t54 = r_i_i_C(1) * t48 - r_i_i_C(2) * t44 + pkin(4);
	t23 = t53 * t41;
	t25 = t53 * t43;
	t52 = t23 * t63 + t25 * t30 - t31 * t35;
	t32 = -t51 * t46 - t50 * t58;
	t33 = -t46 * t58 + t51 * t50;
	t10 = -t24 * t64 - t26 * t32 + t33 * t53;
	t39 = pkin(3) * t49 + pkin(2);
	t29 = -t42 * t50 * t41 + t43 * t61;
	t28 = -t41 * t62 + t43 * t67;
	t22 = -t32 * t41 + t43 * t64;
	t17 = t26 * t33 + t32 * t53;
	t15 = t26 * t31 - t30 * t53;
	t9 = t23 * t64 + t25 * t32 + t33 * t35;
	t2 = t10 * t48 + t22 * t44;
	t1 = -t10 * t44 + t22 * t48;
	t3 = [-t47 * pkin(1) + t7 * pkin(4) + r_i_i_C(1) * t69 + r_i_i_C(2) * t70 + t30 * t28 - t31 * t39 + t51 * t59 + t68 * t52, (t17 * t48 + t33 * t66) * r_i_i_C(1) + (-t17 * t44 + t33 * t65) * r_i_i_C(2) + t17 * pkin(4) + t32 * t39 - t33 * t28 + t68 * (-t25 * t33 + t32 * t35), -t68 * t10 + t54 * t9 + (-t33 * t45 + (t32 * t43 + t41 * t64) * t49) * pkin(3), t22, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t51 * pkin(1) + t10 * pkin(4) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t32 * t28 + t33 * t39 + t47 * t59 + t68 * t9, (t15 * t48 + t31 * t66) * r_i_i_C(1) + (-t15 * t44 + t31 * t65) * r_i_i_C(2) + t15 * pkin(4) - t30 * t39 - t31 * t28 + t68 * (-t25 * t31 - t30 * t35), t68 * t7 - t54 * t52 + (-t31 * t45 + (-t30 * t43 - t41 * t63) * t49) * pkin(3), -t20, -r_i_i_C(1) * t70 + r_i_i_C(2) * t69, 0; 0, (t54 * (t26 * t46 + t50 * t53) + t50 * t39 + (-t28 + (r_i_i_C(1) * t44 + r_i_i_C(2) * t48) * t41) * t46 + t68 * (-t25 * t46 + t35 * t50)) * t42, -t68 * t13 + t54 * (t61 * t23 + (t25 * t50 + t35 * t46) * t42) + (t61 * t41 * t49 + (t43 * t49 * t50 - t45 * t46) * t42) * pkin(3), t29, (-t13 * t44 + t29 * t48) * r_i_i_C(1) + (-t13 * t48 - t29 * t44) * r_i_i_C(2), 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:08
	% DurationCPUTime: 0.90s
	% Computational Cost: add. (842->121), mult. (2275->205), div. (0->0), fcn. (3013->16), ass. (0->72)
	t57 = sin(pkin(7));
	t56 = sin(pkin(13));
	t62 = sin(qJ(3));
	t67 = cos(qJ(3));
	t81 = cos(pkin(13));
	t71 = -t56 * t62 + t67 * t81;
	t39 = t71 * t57;
	t59 = cos(pkin(7));
	t41 = t71 * t59;
	t63 = sin(qJ(2));
	t64 = sin(qJ(1));
	t68 = cos(qJ(2));
	t69 = cos(qJ(1));
	t82 = cos(pkin(6));
	t77 = t69 * t82;
	t46 = t63 * t64 - t68 * t77;
	t47 = t63 * t77 + t64 * t68;
	t51 = -t56 * t67 - t62 * t81;
	t58 = sin(pkin(6));
	t84 = t58 * t69;
	t18 = t39 * t84 + t41 * t46 - t47 * t51;
	t40 = t51 * t57;
	t42 = t51 * t59;
	t19 = -t40 * t84 - t42 * t46 - t47 * t71;
	t35 = -t46 * t57 + t59 * t84;
	t61 = sin(qJ(5));
	t66 = cos(qJ(5));
	t6 = t19 * t66 + t35 * t61;
	t60 = sin(qJ(6));
	t65 = cos(qJ(6));
	t96 = t18 * t65 + t6 * t60;
	t95 = -t18 * t60 + t6 * t65;
	t92 = t19 * t61 - t35 * t66;
	t25 = (-t42 * t68 + t63 * t71) * t58 - t82 * t40;
	t91 = r_i_i_C(3) + pkin(12);
	t90 = pkin(3) * t62;
	t89 = t57 * t61;
	t88 = t57 * t66;
	t87 = t58 * t63;
	t86 = t58 * t64;
	t85 = t58 * t68;
	t83 = pkin(10) + qJ(4);
	t80 = t57 * t87;
	t79 = t58 * (t57 * t90 + t59 * t83 + pkin(9));
	t78 = t64 * t82;
	t74 = r_i_i_C(1) * t65 - r_i_i_C(2) * t60 + pkin(5);
	t73 = -r_i_i_C(1) * t60 - r_i_i_C(2) * t65 - pkin(11);
	t48 = -t63 * t69 - t68 * t78;
	t72 = -t48 * t57 + t59 * t86;
	t49 = -t63 * t78 + t68 * t69;
	t22 = -t40 * t86 - t42 * t48 + t49 * t71;
	t70 = t61 * t91 + t66 * t74 + pkin(4);
	t55 = pkin(3) * t67 + pkin(2);
	t45 = -t57 * t85 + t59 * t82;
	t44 = -t57 * t83 + t59 * t90;
	t33 = (t42 * t63 + t68 * t71) * t58;
	t32 = -t41 * t87 + t51 * t85;
	t31 = t33 * t66 + t61 * t80;
	t29 = t42 * t49 + t48 * t71;
	t28 = -t41 * t49 + t48 * t51;
	t27 = t42 * t47 - t46 * t71;
	t26 = -t41 * t47 - t46 * t51;
	t24 = t82 * t39 + (t41 * t68 + t51 * t63) * t58;
	t21 = t39 * t86 + t41 * t48 + t49 * t51;
	t14 = t29 * t66 + t49 * t89;
	t12 = t27 * t66 + t47 * t89;
	t10 = t25 * t66 + t45 * t61;
	t8 = t22 * t66 + t61 * t72;
	t7 = t22 * t61 - t66 * t72;
	t2 = -t21 * t60 + t65 * t8;
	t1 = -t21 * t65 - t60 * t8;
	t3 = [-t64 * pkin(1) + t19 * pkin(4) + t6 * pkin(5) - t18 * pkin(11) + r_i_i_C(1) * t95 - r_i_i_C(2) * t96 + t46 * t44 - t47 * t55 + t69 * t79 + t91 * t92, (t14 * t65 - t28 * t60) * r_i_i_C(1) + (-t14 * t60 - t28 * t65) * r_i_i_C(2) + t14 * pkin(5) + t29 * pkin(4) - t28 * pkin(11) + t48 * t55 - t49 * t44 + t91 * (t29 * t61 - t49 * t88), -t73 * t22 + (-t49 * t62 + (t48 * t59 + t57 * t86) * t67) * pkin(3) + t70 * t21, t72, -t7 * t74 + t8 * t91, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t69 * pkin(1) + t22 * pkin(4) + t8 * pkin(5) - t21 * pkin(11) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t48 * t44 + t49 * t55 + t64 * t79 + t7 * t91, (t12 * t65 - t26 * t60) * r_i_i_C(1) + (-t12 * t60 - t26 * t65) * r_i_i_C(2) + t12 * pkin(5) + t27 * pkin(4) - t26 * pkin(11) - t46 * t55 - t47 * t44 + t91 * (t27 * t61 - t47 * t88), t73 * t19 + (-t47 * t62 + (-t46 * t59 - t57 * t84) * t67) * pkin(3) - t70 * t18, -t35, -t6 * t91 + t74 * t92, r_i_i_C(1) * t96 + r_i_i_C(2) * t95; 0, (t31 * t65 - t32 * t60) * r_i_i_C(1) + (-t31 * t60 - t32 * t65) * r_i_i_C(2) + t31 * pkin(5) + t33 * pkin(4) - t32 * pkin(11) + (-t44 * t63 + t55 * t68) * t58 + t91 * (t33 * t61 - t66 * t80), -t73 * t25 + (t82 * t57 * t67 + (t59 * t67 * t68 - t62 * t63) * t58) * pkin(3) + t70 * t24, t45, t91 * t10 + t74 * (-t25 * t61 + t45 * t66), (-t10 * t60 - t24 * t65) * r_i_i_C(1) + (-t10 * t65 + t24 * t60) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end