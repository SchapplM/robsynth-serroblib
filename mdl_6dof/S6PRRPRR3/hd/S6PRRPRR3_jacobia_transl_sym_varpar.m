% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:29
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobia_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(6));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(12));
	t1 = sin(pkin(12));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (55->24), mult. (156->50), div. (0->0), fcn. (196->10), ass. (0->24)
	t7 = sin(pkin(7));
	t8 = sin(pkin(6));
	t24 = t7 * t8;
	t10 = cos(pkin(7));
	t15 = cos(qJ(2));
	t23 = t10 * t15;
	t11 = cos(pkin(6));
	t13 = sin(qJ(2));
	t22 = t11 * t13;
	t21 = t11 * t15;
	t12 = sin(qJ(3));
	t14 = cos(qJ(3));
	t20 = r_i_i_C(1) * t14 - r_i_i_C(2) * t12;
	t6 = sin(pkin(12));
	t9 = cos(pkin(12));
	t1 = -t6 * t13 + t9 * t21;
	t19 = -t1 * t10 + t9 * t24;
	t3 = -t9 * t13 - t6 * t21;
	t18 = t10 * t3 + t6 * t24;
	t17 = pkin(2) + t20;
	t16 = (pkin(9) + r_i_i_C(3)) * t7 + (-t12 * r_i_i_C(1) - t14 * r_i_i_C(2)) * t10;
	t4 = t15 * t9 - t6 * t22;
	t2 = t15 * t6 + t9 * t22;
	t5 = [0, t16 * t4 + t17 * t3, (-t12 * t4 + t18 * t14) * r_i_i_C(1) + (-t18 * t12 - t14 * t4) * r_i_i_C(2), 0, 0, 0; 0, t17 * t1 + t16 * t2, (-t12 * t2 - t19 * t14) * r_i_i_C(1) + (t19 * t12 - t14 * t2) * r_i_i_C(2), 0, 0, 0; 1, (t16 * t13 + t17 * t15) * t8, t20 * t7 * t11 + ((-t12 * t13 + t14 * t23) * r_i_i_C(1) + (-t12 * t23 - t13 * t14) * r_i_i_C(2)) * t8, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (112->44), mult. (301->84), div. (0->0), fcn. (383->12), ass. (0->32)
	t25 = cos(qJ(3));
	t36 = t25 * pkin(3);
	t16 = sin(pkin(12));
	t18 = sin(pkin(6));
	t35 = t16 * t18;
	t17 = sin(pkin(7));
	t34 = t17 * t18;
	t20 = cos(pkin(12));
	t33 = t18 * t20;
	t21 = cos(pkin(7));
	t32 = t18 * t21;
	t22 = cos(pkin(6));
	t24 = sin(qJ(2));
	t31 = t22 * t24;
	t26 = cos(qJ(2));
	t30 = t22 * t26;
	t15 = sin(pkin(13));
	t19 = cos(pkin(13));
	t23 = sin(qJ(3));
	t29 = t25 * t15 + t23 * t19;
	t10 = t23 * t15 - t25 * t19;
	t28 = -t10 * r_i_i_C(1) - r_i_i_C(2) * t29 + pkin(2) + t36;
	t3 = t10 * t21;
	t4 = t29 * t21;
	t27 = -t21 * t23 * pkin(3) - t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + (r_i_i_C(3) + pkin(9) + qJ(4)) * t17;
	t9 = -t16 * t31 + t20 * t26;
	t8 = -t16 * t30 - t20 * t24;
	t7 = t16 * t26 + t20 * t31;
	t6 = -t16 * t24 + t20 * t30;
	t2 = t29 * t17;
	t1 = t10 * t17;
	t5 = [0, t27 * t9 + t28 * t8, (-t1 * t35 - t29 * t9 - t8 * t3) * r_i_i_C(1) + (t9 * t10 - t2 * t35 - t8 * t4) * r_i_i_C(2) + (-t9 * t23 + (t16 * t34 + t21 * t8) * t25) * pkin(3), t16 * t32 - t8 * t17, 0, 0; 0, t27 * t7 + t28 * t6, (t1 * t33 - t29 * t7 - t6 * t3) * r_i_i_C(1) + (t7 * t10 + t2 * t33 - t6 * t4) * r_i_i_C(2) + (-t7 * t23 + (-t17 * t33 + t21 * t6) * t25) * pkin(3), -t6 * t17 - t20 * t32, 0, 0; 1, (t27 * t24 + t28 * t26) * t18, (-t1 * r_i_i_C(1) - t2 * r_i_i_C(2) + t17 * t36) * t22 + ((-t24 * t29 - t26 * t3) * r_i_i_C(1) + (t10 * t24 - t26 * t4) * r_i_i_C(2) + (t26 * t21 * t25 - t24 * t23) * pkin(3)) * t18, t22 * t21 - t26 * t34, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (288->73), mult. (781->137), div. (0->0), fcn. (1021->14), ass. (0->45)
	t33 = sin(pkin(13));
	t41 = sin(qJ(3));
	t44 = cos(qJ(3));
	t50 = cos(pkin(13));
	t29 = -t44 * t33 - t41 * t50;
	t35 = sin(pkin(7));
	t19 = t29 * t35;
	t38 = cos(pkin(7));
	t21 = t29 * t38;
	t36 = sin(pkin(6));
	t39 = cos(pkin(6));
	t42 = sin(qJ(2));
	t45 = cos(qJ(2));
	t47 = -t41 * t33 + t44 * t50;
	t9 = -t39 * t19 + (-t21 * t45 + t42 * t47) * t36;
	t60 = -r_i_i_C(3) - pkin(10);
	t34 = sin(pkin(12));
	t59 = t34 * t36;
	t58 = t35 * t36;
	t40 = sin(qJ(5));
	t57 = t35 * t40;
	t43 = cos(qJ(5));
	t56 = t35 * t43;
	t37 = cos(pkin(12));
	t55 = t36 * t37;
	t54 = t36 * t38;
	t52 = t39 * t42;
	t51 = t39 * t45;
	t48 = t43 * r_i_i_C(1) - t40 * r_i_i_C(2) + pkin(4);
	t24 = -t34 * t42 + t37 * t51;
	t25 = t34 * t45 + t37 * t52;
	t46 = -t19 * t55 + t24 * t21 - t25 * t47;
	t26 = -t34 * t51 - t37 * t42;
	t27 = -t34 * t52 + t37 * t45;
	t6 = -t19 * t59 - t26 * t21 + t27 * t47;
	t32 = t44 * pkin(3) + pkin(2);
	t23 = t39 * t38 - t45 * t58;
	t22 = t38 * t41 * pkin(3) + (-pkin(9) - qJ(4)) * t35;
	t20 = t47 * t38;
	t18 = t47 * t35;
	t17 = -t26 * t35 + t34 * t54;
	t16 = -t24 * t35 - t37 * t54;
	t13 = t27 * t21 + t26 * t47;
	t11 = t25 * t21 + t24 * t47;
	t1 = [0, (t13 * t43 + t27 * t57) * r_i_i_C(1) + (-t13 * t40 + t27 * t56) * r_i_i_C(2) + t13 * pkin(4) + t26 * t32 - t27 * t22 + t60 * (-t27 * t20 + t26 * t29), -t60 * t6 + t48 * (t18 * t59 + t26 * t20 + t27 * t29) + (-t27 * t41 + (t26 * t38 + t34 * t58) * t44) * pkin(3), t17, (t17 * t43 - t6 * t40) * r_i_i_C(1) + (-t17 * t40 - t6 * t43) * r_i_i_C(2), 0; 0, (t11 * t43 + t25 * t57) * r_i_i_C(1) + (-t11 * t40 + t25 * t56) * r_i_i_C(2) + t11 * pkin(4) + t24 * t32 - t25 * t22 + t60 * (-t25 * t20 + t24 * t29), t60 * t46 + t48 * (-t18 * t55 + t24 * t20 + t25 * t29) + (-t25 * t41 + (t24 * t38 - t35 * t55) * t44) * pkin(3), t16, (t16 * t43 + t40 * t46) * r_i_i_C(1) + (-t16 * t40 + t43 * t46) * r_i_i_C(2), 0; 1, (t48 * (t21 * t42 + t45 * t47) + t45 * t32 + (-t22 + (t40 * r_i_i_C(1) + t43 * r_i_i_C(2)) * t35) * t42 + t60 * (-t20 * t42 + t29 * t45)) * t36, -t60 * t9 + t48 * (t39 * t18 + (t20 * t45 + t29 * t42) * t36) + (t35 * t39 * t44 + (t38 * t44 * t45 - t41 * t42) * t36) * pkin(3), t23, (t23 * t43 - t9 * t40) * r_i_i_C(1) + (-t23 * t40 - t9 * t43) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.69s
	% Computational Cost: add. (646->100), mult. (1769->183), div. (0->0), fcn. (2343->16), ass. (0->64)
	t48 = sin(pkin(13));
	t57 = sin(qJ(3));
	t61 = cos(qJ(3));
	t70 = cos(pkin(13));
	t44 = -t61 * t48 - t57 * t70;
	t50 = sin(pkin(7));
	t34 = t44 * t50;
	t53 = cos(pkin(7));
	t36 = t44 * t53;
	t51 = sin(pkin(6));
	t54 = cos(pkin(6));
	t58 = sin(qJ(2));
	t62 = cos(qJ(2));
	t65 = -t57 * t48 + t61 * t70;
	t19 = -t54 * t34 + (-t36 * t62 + t58 * t65) * t51;
	t81 = r_i_i_C(3) + pkin(11);
	t49 = sin(pkin(12));
	t80 = t49 * t51;
	t56 = sin(qJ(5));
	t79 = t50 * t56;
	t60 = cos(qJ(5));
	t78 = t50 * t60;
	t52 = cos(pkin(12));
	t77 = t51 * t52;
	t76 = t51 * t53;
	t75 = t51 * t58;
	t74 = t51 * t62;
	t72 = t54 * t58;
	t71 = t54 * t62;
	t69 = t50 * t75;
	t55 = sin(qJ(6));
	t59 = cos(qJ(6));
	t67 = t59 * r_i_i_C(1) - t55 * r_i_i_C(2) + pkin(5);
	t66 = -t55 * r_i_i_C(1) - t59 * r_i_i_C(2) - pkin(10);
	t39 = -t49 * t58 + t52 * t71;
	t40 = t49 * t62 + t52 * t72;
	t64 = -t34 * t77 + t39 * t36 - t40 * t65;
	t41 = -t49 * t71 - t52 * t58;
	t42 = -t49 * t72 + t52 * t62;
	t14 = -t34 * t80 - t41 * t36 + t42 * t65;
	t63 = t81 * t56 + t67 * t60 + pkin(4);
	t47 = t61 * pkin(3) + pkin(2);
	t38 = -t50 * t74 + t54 * t53;
	t37 = t53 * t57 * pkin(3) + (-pkin(9) - qJ(4)) * t50;
	t35 = t65 * t53;
	t33 = t65 * t50;
	t29 = -t41 * t50 + t49 * t76;
	t28 = -t39 * t50 - t52 * t76;
	t27 = (t36 * t58 + t62 * t65) * t51;
	t26 = -t35 * t75 + t44 * t74;
	t25 = t27 * t60 + t56 * t69;
	t23 = t42 * t36 + t41 * t65;
	t22 = -t42 * t35 + t41 * t44;
	t21 = t40 * t36 + t39 * t65;
	t20 = -t40 * t35 + t39 * t44;
	t18 = t54 * t33 + (t35 * t62 + t44 * t58) * t51;
	t16 = t19 * t60 + t38 * t56;
	t13 = t33 * t80 + t41 * t35 + t42 * t44;
	t10 = -t33 * t77 + t39 * t35 + t40 * t44;
	t8 = t23 * t60 + t42 * t79;
	t6 = t21 * t60 + t40 * t79;
	t4 = t14 * t60 + t29 * t56;
	t2 = t28 * t56 - t60 * t64;
	t1 = [0, (-t22 * t55 + t8 * t59) * r_i_i_C(1) + (-t22 * t59 - t8 * t55) * r_i_i_C(2) + t8 * pkin(5) + t23 * pkin(4) - t22 * pkin(10) + t41 * t47 - t42 * t37 + t81 * (t23 * t56 - t42 * t78), -t66 * t14 + (-t42 * t57 + (t41 * t53 + t50 * t80) * t61) * pkin(3) + t63 * t13, t29, t81 * t4 + t67 * (-t14 * t56 + t29 * t60), (-t13 * t59 - t4 * t55) * r_i_i_C(1) + (t13 * t55 - t4 * t59) * r_i_i_C(2); 0, (-t20 * t55 + t6 * t59) * r_i_i_C(1) + (-t20 * t59 - t6 * t55) * r_i_i_C(2) + t6 * pkin(5) + t21 * pkin(4) - t20 * pkin(10) + t39 * t47 - t40 * t37 + t81 * (t21 * t56 - t40 * t78), t66 * t64 + (-t40 * t57 + (t39 * t53 - t50 * t77) * t61) * pkin(3) + t63 * t10, t28, t81 * t2 + t67 * (t28 * t60 + t56 * t64), (-t10 * t59 - t2 * t55) * r_i_i_C(1) + (t10 * t55 - t2 * t59) * r_i_i_C(2); 1, (t25 * t59 - t26 * t55) * r_i_i_C(1) + (-t25 * t55 - t26 * t59) * r_i_i_C(2) + t25 * pkin(5) + t27 * pkin(4) - t26 * pkin(10) + (-t58 * t37 + t62 * t47) * t51 + t81 * (t27 * t56 - t60 * t69), -t66 * t19 + (t54 * t50 * t61 + (t53 * t61 * t62 - t57 * t58) * t51) * pkin(3) + t63 * t18, t38, t81 * t16 + t67 * (-t19 * t56 + t38 * t60), (-t16 * t55 - t18 * t59) * r_i_i_C(1) + (-t16 * t59 + t18 * t55) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end