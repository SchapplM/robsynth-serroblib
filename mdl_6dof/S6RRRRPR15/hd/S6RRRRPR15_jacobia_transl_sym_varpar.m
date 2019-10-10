% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:56
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR15_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR15_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR15_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
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
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
	% DurationCPUTime: 0.13s
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
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
	% DurationCPUTime: 0.24s
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
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
	% DurationCPUTime: 0.42s
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
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (390->78), mult. (1067->135), div. (0->0), fcn. (1387->12), ass. (0->55)
	t47 = sin(qJ(2));
	t50 = cos(qJ(2));
	t51 = cos(qJ(1));
	t62 = cos(pkin(6));
	t58 = t51 * t62;
	t74 = sin(qJ(1));
	t34 = t74 * t47 - t50 * t58;
	t35 = t47 * t58 + t74 * t50;
	t46 = sin(qJ(3));
	t49 = cos(qJ(3));
	t42 = sin(pkin(7));
	t43 = sin(pkin(6));
	t70 = t43 * t51;
	t60 = t42 * t70;
	t44 = cos(pkin(7));
	t69 = t44 * t46;
	t16 = t34 * t69 - t35 * t49 + t46 * t60;
	t28 = -t34 * t42 + t44 * t70;
	t45 = sin(qJ(4));
	t48 = cos(qJ(4));
	t79 = t16 * t48 + t28 * t45;
	t78 = t16 * t45 - t28 * t48;
	t63 = r_i_i_C(3) + qJ(5);
	t77 = pkin(4) - r_i_i_C(2);
	t53 = t63 * t45 + t77 * t48 + pkin(3);
	t76 = pkin(11) + r_i_i_C(1);
	t75 = pkin(10) * t42;
	t73 = t42 * t45;
	t72 = t42 * t47;
	t71 = t42 * t48;
	t68 = t44 * t49;
	t67 = t46 * t47;
	t66 = t46 * t50;
	t65 = t47 * t49;
	t64 = t49 * t50;
	t61 = t43 * t72;
	t59 = t43 * t74;
	t57 = t62 * t42;
	t56 = t42 * t59;
	t55 = t62 * t74;
	t36 = -t51 * t47 - t50 * t55;
	t54 = -t36 * t42 + t44 * t59;
	t52 = -t35 * t46 + (-t34 * t44 - t60) * t49;
	t37 = -t47 * t55 + t51 * t50;
	t33 = -t43 * t50 * t42 + t62 * t44;
	t32 = (-t44 * t67 + t64) * t43;
	t27 = t46 * t57 + (t44 * t66 + t65) * t43;
	t22 = t36 * t49 - t37 * t69;
	t20 = -t34 * t49 - t35 * t69;
	t18 = t37 * t49 + (t36 * t44 + t56) * t46;
	t17 = -t36 * t68 + t37 * t46 - t49 * t56;
	t11 = t27 * t45 - t33 * t48;
	t6 = t18 * t48 + t54 * t45;
	t5 = t18 * t45 - t54 * t48;
	t1 = [-t74 * pkin(1) - t35 * pkin(2) + t16 * pkin(3) + pkin(9) * t70 + t28 * pkin(10) + t76 * t52 + t63 * t78 + t77 * t79, t37 * t75 + t36 * pkin(2) + t22 * pkin(3) + t63 * (t22 * t45 - t37 * t71) + t76 * (t36 * t46 + t37 * t68) + t77 * (t22 * t48 + t37 * t73), -t53 * t17 + t76 * t18, -t77 * t5 + t63 * t6, t5, 0; t51 * pkin(1) + t37 * pkin(2) + t18 * pkin(3) + pkin(9) * t59 + t54 * pkin(10) + t76 * t17 + t63 * t5 + t77 * t6, t35 * t75 - t34 * pkin(2) + t20 * pkin(3) + t77 * (t20 * t48 + t35 * t73) + t63 * (t20 * t45 - t35 * t71) + t76 * (-t34 * t46 + t35 * t68), -t16 * t76 + t53 * t52, -t63 * t79 + t77 * t78, -t78, 0; 0, t32 * pkin(3) + t77 * (t32 * t48 + t45 * t61) + t63 * (t32 * t45 - t48 * t61) + (pkin(2) * t50 + pkin(10) * t72 + t76 * (t44 * t65 + t66)) * t43, t76 * t27 + t53 * (t49 * t57 + (t44 * t64 - t67) * t43), t63 * (t27 * t48 + t33 * t45) - t77 * t11, t11, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (651->90), mult. (1797->159), div. (0->0), fcn. (2353->14), ass. (0->59)
	t50 = sin(qJ(2));
	t53 = cos(qJ(2));
	t54 = cos(qJ(1));
	t75 = cos(pkin(6));
	t67 = t54 * t75;
	t81 = sin(qJ(1));
	t37 = t81 * t50 - t53 * t67;
	t38 = t50 * t67 + t81 * t53;
	t49 = sin(qJ(3));
	t74 = cos(pkin(7));
	t68 = t49 * t74;
	t45 = sin(pkin(7));
	t46 = sin(pkin(6));
	t77 = t46 * t54;
	t71 = t45 * t77;
	t82 = cos(qJ(3));
	t16 = -t37 * t68 + t38 * t82 - t49 * t71;
	t48 = sin(qJ(4));
	t52 = cos(qJ(4));
	t69 = t46 * t74;
	t60 = t37 * t45 - t54 * t69;
	t86 = t16 * t52 + t60 * t48;
	t3 = t16 * t48 - t60 * t52;
	t64 = t75 * t81;
	t59 = t54 * t50 + t53 * t64;
	t70 = t46 * t81;
	t85 = -t45 * t70 + t59 * t74;
	t84 = pkin(5) + pkin(11);
	t83 = pkin(10) * t45;
	t80 = t45 * t48;
	t79 = t45 * t52;
	t78 = t46 * t53;
	t76 = t50 * t45;
	t73 = r_i_i_C(3) + pkin(12) + pkin(4);
	t72 = t46 * t76;
	t66 = t75 * t45;
	t63 = t74 * t82;
	t47 = sin(qJ(6));
	t51 = cos(qJ(6));
	t62 = t47 * r_i_i_C(1) + t51 * r_i_i_C(2) + qJ(5);
	t61 = t51 * r_i_i_C(1) - t47 * r_i_i_C(2) + t84;
	t58 = -t45 * t78 + t75 * t74;
	t56 = -t62 * t48 - t73 * t52 - pkin(3);
	t15 = t37 * t63 + t38 * t49 + t82 * t71;
	t55 = t59 * t45 + t81 * t69;
	t39 = -t50 * t64 + t54 * t53;
	t35 = (-t50 * t68 + t82 * t53) * t46;
	t30 = t49 * t66 + (t82 * t50 + t53 * t68) * t46;
	t29 = t46 * t50 * t49 - t63 * t78 - t82 * t66;
	t24 = -t39 * t68 - t59 * t82;
	t22 = -t37 * t82 - t38 * t68;
	t20 = t39 * t82 - t85 * t49;
	t19 = t39 * t49 + t85 * t82;
	t13 = t30 * t48 - t58 * t52;
	t8 = t20 * t52 + t55 * t48;
	t7 = t20 * t48 - t55 * t52;
	t2 = t19 * t51 + t7 * t47;
	t1 = -t19 * t47 + t7 * t51;
	t4 = [-t81 * pkin(1) - t38 * pkin(2) - t16 * pkin(3) + pkin(9) * t77 - t60 * pkin(10) - t61 * t15 - t62 * t3 - t73 * t86, t24 * pkin(3) - t59 * pkin(2) + t39 * t83 + t62 * (t24 * t48 - t39 * t79) + t61 * (t39 * t63 - t59 * t49) + t73 * (t24 * t52 + t39 * t80), t56 * t19 + t61 * t20, t62 * t8 - t73 * t7, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t54 * pkin(1) + t39 * pkin(2) + t20 * pkin(3) + pkin(9) * t70 + t55 * pkin(10) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(5) + t84 * t19 + t73 * t8, t38 * t83 - t37 * pkin(2) + t22 * pkin(3) + t62 * (t22 * t48 - t38 * t79) + t61 * (-t37 * t49 + t38 * t63) + t73 * (t22 * t52 + t38 * t80), t56 * t15 + t61 * t16, -t73 * t3 + t62 * t86, t3, (-t15 * t47 + t3 * t51) * r_i_i_C(1) + (-t15 * t51 - t3 * t47) * r_i_i_C(2); 0, t35 * pkin(3) + t62 * (t35 * t48 - t52 * t72) + t73 * (t35 * t52 + t48 * t72) + (t53 * pkin(2) + pkin(10) * t76 + t61 * (t49 * t53 + t50 * t63)) * t46, t56 * t29 + t61 * t30, t62 * (t30 * t52 + t58 * t48) - t73 * t13, t13, (t13 * t51 - t29 * t47) * r_i_i_C(1) + (-t13 * t47 - t29 * t51) * r_i_i_C(2);];
	Ja_transl = t4;
else
	Ja_transl=NaN(3,6);
end