% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRR15
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR15_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR15_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR15_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
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
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
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
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
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
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (169->48), mult. (455->84), div. (0->0), fcn. (582->10), ass. (0->37)
	t54 = pkin(10) + r_i_i_C(1);
	t27 = sin(qJ(2));
	t28 = sin(qJ(1));
	t30 = cos(qJ(2));
	t31 = cos(qJ(1));
	t40 = cos(pkin(6));
	t36 = t31 * t40;
	t17 = t27 * t36 + t28 * t30;
	t26 = sin(qJ(3));
	t29 = cos(qJ(3));
	t16 = t28 * t27 - t30 * t36;
	t23 = sin(pkin(7));
	t25 = cos(pkin(7));
	t24 = sin(pkin(6));
	t48 = t24 * t31;
	t33 = t16 * t25 + t23 * t48;
	t53 = -t17 * t29 + t26 * t33;
	t1 = t17 * t26 + t29 * t33;
	t52 = -r_i_i_C(2) + pkin(3);
	t49 = t24 * t28;
	t47 = t25 * t26;
	t46 = t25 * t29;
	t45 = t26 * t27;
	t44 = t26 * t30;
	t43 = t27 * t29;
	t42 = t29 * t30;
	t41 = r_i_i_C(3) + qJ(4);
	t39 = t23 * t49;
	t38 = t54 * t23;
	t37 = t28 * t40;
	t35 = t40 * t23;
	t19 = -t27 * t37 + t31 * t30;
	t18 = -t31 * t27 - t30 * t37;
	t11 = -t29 * t35 + (-t25 * t42 + t45) * t24;
	t6 = t19 * t29 + (t18 * t25 + t39) * t26;
	t5 = -t18 * t46 + t19 * t26 - t29 * t39;
	t2 = [-t17 * pkin(2) - t28 * pkin(1) + pkin(9) * t48 + t52 * t53 - t41 * t1 + t54 * (-t16 * t23 + t25 * t48), t18 * pkin(2) + t41 * (t18 * t26 + t19 * t46) + t19 * t38 + t52 * (t18 * t29 - t19 * t47), t41 * t6 - t52 * t5, t5, 0, 0; t31 * pkin(1) + t19 * pkin(2) + pkin(9) * t49 + t41 * t5 + t52 * t6 + t54 * (-t18 * t23 + t25 * t49), -t16 * pkin(2) + t52 * (-t16 * t29 - t17 * t47) + t41 * (-t16 * t26 + t17 * t46) + t17 * t38, -t52 * t1 - t41 * t53, t1, 0, 0; 0, (t52 * (-t25 * t45 + t42) + t41 * (t25 * t43 + t44) + pkin(2) * t30 + t27 * t38) * t24, t41 * (t26 * t35 + (t25 * t44 + t43) * t24) - t52 * t11, t11, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (296->65), mult. (807->113), div. (0->0), fcn. (1042->12), ass. (0->47)
	t65 = pkin(4) + pkin(10);
	t37 = sin(qJ(2));
	t38 = sin(qJ(1));
	t41 = cos(qJ(2));
	t42 = cos(qJ(1));
	t52 = cos(pkin(6));
	t47 = t42 * t52;
	t24 = t37 * t47 + t38 * t41;
	t36 = sin(qJ(3));
	t40 = cos(qJ(3));
	t23 = t38 * t37 - t41 * t47;
	t34 = cos(pkin(7));
	t32 = sin(pkin(7));
	t33 = sin(pkin(6));
	t59 = t33 * t42;
	t49 = t32 * t59;
	t44 = t23 * t34 + t49;
	t64 = -t24 * t40 + t44 * t36;
	t35 = sin(qJ(5));
	t39 = cos(qJ(5));
	t45 = t35 * r_i_i_C(1) + t39 * r_i_i_C(2) + qJ(4);
	t63 = t39 * r_i_i_C(1) - t35 * r_i_i_C(2) + t65;
	t62 = t24 * t36;
	t60 = t33 * t38;
	t58 = t34 * t36;
	t57 = t34 * t40;
	t56 = t36 * t37;
	t55 = t36 * t41;
	t54 = t37 * t40;
	t53 = t40 * t41;
	t51 = r_i_i_C(3) + pkin(11) + pkin(3);
	t50 = t32 * t60;
	t48 = t38 * t52;
	t46 = t52 * t32;
	t15 = -t23 * t32 + t34 * t59;
	t25 = -t42 * t37 - t41 * t48;
	t17 = -t25 * t32 + t34 * t60;
	t43 = t63 * t32;
	t26 = -t37 * t48 + t42 * t41;
	t22 = -t33 * t41 * t32 + t52 * t34;
	t13 = -t40 * t46 + (-t34 * t53 + t56) * t33;
	t8 = t26 * t40 + (t25 * t34 + t50) * t36;
	t7 = -t25 * t57 + t26 * t36 - t40 * t50;
	t3 = t23 * t57 + t40 * t49 + t62;
	t2 = t17 * t39 + t7 * t35;
	t1 = -t17 * t35 + t7 * t39;
	t4 = [-t24 * pkin(2) - t38 * pkin(1) + pkin(9) * t59 + t51 * t64 + t45 * (-t44 * t40 - t62) + t63 * t15, t25 * pkin(2) + t45 * (t25 * t36 + t26 * t57) + t26 * t43 + t51 * (t25 * t40 - t26 * t58), t45 * t8 - t51 * t7, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t42 * pkin(1) + t26 * pkin(2) + pkin(9) * t60 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(4) + t65 * t17 + t51 * t8, -t23 * pkin(2) + t45 * (-t23 * t36 + t24 * t57) + t24 * t43 + t51 * (-t23 * t40 - t24 * t58), -t51 * t3 - t45 * t64, t3, (t15 * t35 + t3 * t39) * r_i_i_C(1) + (t15 * t39 - t3 * t35) * r_i_i_C(2), 0; 0, (t45 * (t34 * t54 + t55) + t41 * pkin(2) + t37 * t43 + t51 * (-t34 * t56 + t53)) * t33, t45 * (t36 * t46 + (t34 * t55 + t54) * t33) - t51 * t13, t13, (t13 * t39 - t22 * t35) * r_i_i_C(1) + (-t13 * t35 - t22 * t39) * r_i_i_C(2), 0;];
	Ja_transl = t4;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.57s
	% Computational Cost: add. (589->93), mult. (1627->159), div. (0->0), fcn. (2127->14), ass. (0->63)
	t88 = pkin(4) + pkin(10);
	t53 = sin(qJ(2));
	t54 = sin(qJ(1));
	t58 = cos(qJ(2));
	t59 = cos(qJ(1));
	t71 = cos(pkin(6));
	t65 = t59 * t71;
	t39 = t53 * t65 + t54 * t58;
	t52 = sin(qJ(3));
	t57 = cos(qJ(3));
	t38 = t54 * t53 - t58 * t65;
	t49 = cos(pkin(7));
	t47 = sin(pkin(7));
	t48 = sin(pkin(6));
	t78 = t48 * t59;
	t68 = t47 * t78;
	t62 = t38 * t49 + t68;
	t87 = -t39 * t57 + t62 * t52;
	t51 = sin(qJ(5));
	t56 = cos(qJ(5));
	t50 = sin(qJ(6));
	t55 = cos(qJ(6));
	t63 = t55 * r_i_i_C(1) - t50 * r_i_i_C(2) + pkin(5);
	t85 = r_i_i_C(3) + pkin(12);
	t60 = t63 * t51 - t85 * t56 + qJ(4);
	t86 = pkin(3) + pkin(11);
	t84 = t39 * t52;
	t82 = t47 * t48;
	t81 = t47 * t51;
	t80 = t47 * t56;
	t79 = t48 * t54;
	t77 = t49 * t52;
	t76 = t49 * t57;
	t75 = t52 * t53;
	t74 = t52 * t58;
	t73 = t53 * t57;
	t72 = t57 * t58;
	t70 = t53 * t82;
	t69 = t47 * t79;
	t67 = t88 * t47;
	t66 = t54 * t71;
	t64 = t71 * t47;
	t30 = -t38 * t47 + t49 * t78;
	t40 = -t59 * t53 - t58 * t66;
	t32 = -t40 * t47 + t49 * t79;
	t61 = t50 * r_i_i_C(1) + t55 * r_i_i_C(2) + t86;
	t41 = -t53 * t66 + t59 * t58;
	t37 = t71 * t49 - t58 * t82;
	t35 = (t49 * t73 + t74) * t48;
	t29 = t52 * t64 + (t49 * t74 + t73) * t48;
	t28 = -t57 * t64 + (-t49 * t72 + t75) * t48;
	t24 = t40 * t52 + t41 * t76;
	t22 = -t38 * t52 + t39 * t76;
	t21 = t41 * t57 + (t40 * t49 + t69) * t52;
	t20 = -t40 * t76 + t41 * t52 - t57 * t69;
	t16 = t38 * t76 + t57 * t68 + t84;
	t15 = t28 * t51 + t37 * t56;
	t8 = t20 * t51 + t32 * t56;
	t7 = -t20 * t56 + t32 * t51;
	t6 = t16 * t51 - t30 * t56;
	t2 = t21 * t50 + t8 * t55;
	t1 = t21 * t55 - t8 * t50;
	t3 = [pkin(9) * t78 - t54 * pkin(1) - t39 * pkin(2) + t61 * t87 + t60 * (-t62 * t57 - t84) + (t85 * t51 + t63 * t56 + t88) * t30, t40 * pkin(2) + t24 * qJ(4) + t41 * t67 + t63 * (t24 * t51 + t41 * t80) + t61 * (t40 * t57 - t41 * t77) - t85 * (t24 * t56 - t41 * t81), -t61 * t20 + t60 * t21, t20, -t63 * t7 + t85 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t59 * pkin(1) + t41 * pkin(2) + t8 * pkin(5) + pkin(9) * t79 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t20 * qJ(4) + t86 * t21 + t88 * t32 + t85 * t7, -t38 * pkin(2) + t22 * qJ(4) - t85 * (t22 * t56 - t39 * t81) + t39 * t67 + t63 * (t22 * t51 + t39 * t80) + t61 * (-t38 * t57 - t39 * t77), -t61 * t16 - t60 * t87, t16, t85 * t6 + t63 * (t16 * t56 + t30 * t51), (-t6 * t50 - t55 * t87) * r_i_i_C(1) + (t50 * t87 - t6 * t55) * r_i_i_C(2); 0, t35 * qJ(4) + t63 * (t35 * t51 + t56 * t70) + t85 * (-t35 * t56 + t51 * t70) + (t61 * (-t49 * t75 + t72) + t58 * pkin(2) + t53 * t67) * t48, -t61 * t28 + t60 * t29, t28, t85 * t15 + t63 * (t28 * t56 - t37 * t51), (-t15 * t50 + t29 * t55) * r_i_i_C(1) + (-t15 * t55 - t29 * t50) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end