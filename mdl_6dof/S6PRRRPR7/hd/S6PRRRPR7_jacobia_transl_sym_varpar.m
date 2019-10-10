% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:58
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR7_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR7_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR7_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobia_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.14s
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
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
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
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (168->57), mult. (471->110), div. (0->0), fcn. (606->12), ass. (0->43)
	t50 = r_i_i_C(3) + pkin(10);
	t22 = sin(pkin(7));
	t49 = t22 * pkin(9);
	t23 = sin(pkin(6));
	t48 = t22 * t23;
	t26 = cos(pkin(6));
	t47 = t22 * t26;
	t27 = sin(qJ(4));
	t46 = t22 * t27;
	t30 = cos(qJ(4));
	t45 = t22 * t30;
	t25 = cos(pkin(7));
	t44 = t23 * t25;
	t28 = sin(qJ(3));
	t43 = t25 * t28;
	t31 = cos(qJ(3));
	t42 = t25 * t31;
	t29 = sin(qJ(2));
	t41 = t26 * t29;
	t32 = cos(qJ(2));
	t40 = t26 * t32;
	t39 = t28 * t29;
	t38 = t28 * t32;
	t37 = t29 * t31;
	t36 = t31 * t32;
	t35 = t30 * r_i_i_C(1) - t27 * r_i_i_C(2) + pkin(3);
	t21 = sin(pkin(12));
	t24 = cos(pkin(12));
	t16 = -t21 * t29 + t24 * t40;
	t34 = t16 * t25 - t24 * t48;
	t18 = -t21 * t40 - t24 * t29;
	t33 = t18 * t25 + t21 * t48;
	t19 = -t21 * t41 + t24 * t32;
	t17 = t21 * t32 + t24 * t41;
	t15 = t26 * t25 - t32 * t48;
	t12 = -t18 * t22 + t21 * t44;
	t11 = -t16 * t22 - t24 * t44;
	t10 = t28 * t47 + (t25 * t38 + t37) * t23;
	t8 = t18 * t31 - t19 * t43;
	t6 = t16 * t31 - t17 * t43;
	t4 = t19 * t31 + t33 * t28;
	t2 = t17 * t31 + t34 * t28;
	t1 = [0, (t19 * t46 + t8 * t30) * r_i_i_C(1) + (t19 * t45 - t8 * t27) * r_i_i_C(2) + t8 * pkin(3) + t18 * pkin(2) + t19 * t49 + t50 * (t18 * t28 + t19 * t42), t50 * t4 + t35 * (-t19 * t28 + t33 * t31), (t12 * t30 - t4 * t27) * r_i_i_C(1) + (-t12 * t27 - t4 * t30) * r_i_i_C(2), 0, 0; 0, (t17 * t46 + t6 * t30) * r_i_i_C(1) + (t17 * t45 - t6 * t27) * r_i_i_C(2) + t6 * pkin(3) + t16 * pkin(2) + t17 * t49 + t50 * (t16 * t28 + t17 * t42), t50 * t2 + t35 * (-t17 * t28 + t34 * t31), (t11 * t30 - t2 * t27) * r_i_i_C(1) + (-t11 * t27 - t2 * t30) * r_i_i_C(2), 0, 0; 1, (t35 * (-t25 * t39 + t36) + t32 * pkin(2) + (t27 * r_i_i_C(1) + t30 * r_i_i_C(2) + pkin(9)) * t29 * t22 + t50 * (t25 * t37 + t38)) * t23, t50 * t10 + t35 * (t31 * t47 + (t25 * t36 - t39) * t23), (-t10 * t27 + t15 * t30) * r_i_i_C(1) + (-t10 * t30 - t15 * t27) * r_i_i_C(2), 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (356->75), mult. (1003->139), div. (0->0), fcn. (1308->14), ass. (0->59)
	t35 = sin(pkin(7));
	t71 = pkin(9) * t35;
	t36 = sin(pkin(6));
	t70 = t35 * t36;
	t40 = cos(pkin(6));
	t69 = t35 * t40;
	t41 = sin(qJ(4));
	t68 = t35 * t41;
	t43 = sin(qJ(2));
	t67 = t35 * t43;
	t44 = cos(qJ(4));
	t66 = t35 * t44;
	t39 = cos(pkin(7));
	t65 = t36 * t39;
	t42 = sin(qJ(3));
	t64 = t39 * t42;
	t45 = cos(qJ(3));
	t63 = t39 * t45;
	t62 = t40 * t43;
	t46 = cos(qJ(2));
	t61 = t40 * t46;
	t60 = t42 * t43;
	t59 = t42 * t46;
	t58 = t43 * t45;
	t57 = t45 * t46;
	t56 = r_i_i_C(3) + qJ(5);
	t55 = t36 * t67;
	t33 = sin(pkin(13));
	t37 = cos(pkin(13));
	t54 = r_i_i_C(1) * t37 - r_i_i_C(2) * t33 + pkin(4);
	t53 = t33 * r_i_i_C(1) + t37 * r_i_i_C(2) + pkin(10);
	t34 = sin(pkin(12));
	t38 = cos(pkin(12));
	t28 = -t34 * t43 + t38 * t61;
	t52 = -t28 * t35 - t38 * t65;
	t51 = t28 * t39 - t38 * t70;
	t30 = -t34 * t61 - t38 * t43;
	t50 = -t30 * t35 + t34 * t65;
	t49 = t30 * t39 + t34 * t70;
	t48 = t40 * t39 - t46 * t70;
	t47 = t56 * t41 + t54 * t44 + pkin(3);
	t31 = -t34 * t62 + t38 * t46;
	t29 = t34 * t46 + t38 * t62;
	t26 = (-t39 * t60 + t57) * t36;
	t25 = (t39 * t58 + t59) * t36;
	t24 = t42 * t69 + (t39 * t59 + t58) * t36;
	t20 = t26 * t44 + t41 * t55;
	t18 = t30 * t45 - t31 * t64;
	t17 = t30 * t42 + t31 * t63;
	t16 = t28 * t45 - t29 * t64;
	t15 = t28 * t42 + t29 * t63;
	t13 = t24 * t41 - t48 * t44;
	t12 = t31 * t45 + t49 * t42;
	t10 = t29 * t45 + t51 * t42;
	t8 = t18 * t44 + t31 * t68;
	t6 = t16 * t44 + t29 * t68;
	t3 = t12 * t41 - t50 * t44;
	t1 = t10 * t41 - t52 * t44;
	t2 = [0, (t17 * t33 + t8 * t37) * r_i_i_C(1) + (t17 * t37 - t8 * t33) * r_i_i_C(2) + t8 * pkin(4) + t18 * pkin(3) + t17 * pkin(10) + t30 * pkin(2) + t31 * t71 + t56 * (t18 * t41 - t31 * t66), t53 * t12 + t47 * (-t31 * t42 + t49 * t45), t56 * (t12 * t44 + t50 * t41) - t54 * t3, t3, 0; 0, (t15 * t33 + t6 * t37) * r_i_i_C(1) + (t15 * t37 - t6 * t33) * r_i_i_C(2) + t6 * pkin(4) + t16 * pkin(3) + t15 * pkin(10) + t28 * pkin(2) + t29 * t71 + t56 * (t16 * t41 - t29 * t66), t53 * t10 + t47 * (-t29 * t42 + t51 * t45), t56 * (t10 * t44 + t52 * t41) - t54 * t1, t1, 0; 1, (t20 * t37 + t25 * t33) * r_i_i_C(1) + (-t20 * t33 + t25 * t37) * r_i_i_C(2) + t20 * pkin(4) + t26 * pkin(3) + t25 * pkin(10) + (t46 * pkin(2) + pkin(9) * t67) * t36 + t56 * (t26 * t41 - t44 * t55), t53 * t24 + t47 * (t45 * t69 + (t39 * t57 - t60) * t36), t56 * (t24 * t44 + t48 * t41) - t54 * t13, t13, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (504->75), mult. (1262->140), div. (0->0), fcn. (1646->16), ass. (0->57)
	t41 = sin(pkin(7));
	t47 = sin(qJ(2));
	t49 = cos(qJ(2));
	t43 = cos(pkin(12));
	t71 = cos(pkin(6));
	t65 = t43 * t71;
	t69 = sin(pkin(12));
	t56 = t69 * t47 - t49 * t65;
	t70 = cos(pkin(7));
	t42 = sin(pkin(6));
	t74 = t42 * t43;
	t81 = t41 * t74 + t56 * t70;
	t60 = t71 * t69;
	t57 = t43 * t47 + t49 * t60;
	t66 = t42 * t69;
	t80 = -t41 * t66 + t57 * t70;
	t79 = pkin(9) * t41;
	t78 = r_i_i_C(3) + pkin(11) + qJ(5);
	t77 = cos(qJ(3));
	t45 = sin(qJ(4));
	t76 = t41 * t45;
	t48 = cos(qJ(4));
	t75 = t41 * t48;
	t73 = t42 * t47;
	t72 = t42 * t49;
	t68 = t41 * t73;
	t46 = sin(qJ(3));
	t64 = t46 * t70;
	t63 = t71 * t41;
	t61 = t70 * t77;
	t39 = pkin(13) + qJ(6);
	t37 = sin(t39);
	t38 = cos(t39);
	t59 = t38 * r_i_i_C(1) - t37 * r_i_i_C(2) + cos(pkin(13)) * pkin(5) + pkin(4);
	t58 = sin(pkin(13)) * pkin(5) + t37 * r_i_i_C(1) + t38 * r_i_i_C(2) + pkin(10);
	t55 = -t41 * t72 + t71 * t70;
	t54 = -t78 * t45 - t59 * t48 - pkin(3);
	t51 = t56 * t41 - t70 * t74;
	t50 = t57 * t41 + t70 * t66;
	t31 = t43 * t49 - t47 * t60;
	t30 = t47 * t65 + t69 * t49;
	t28 = (-t47 * t64 + t77 * t49) * t42;
	t24 = t46 * t63 + (t77 * t47 + t49 * t64) * t42;
	t23 = t46 * t73 - t61 * t72 - t77 * t63;
	t18 = -t31 * t64 - t57 * t77;
	t16 = -t30 * t64 - t56 * t77;
	t14 = t24 * t48 + t55 * t45;
	t13 = t24 * t45 - t55 * t48;
	t12 = t31 * t77 - t80 * t46;
	t11 = t31 * t46 + t80 * t77;
	t10 = t30 * t77 - t81 * t46;
	t9 = t30 * t46 + t81 * t77;
	t4 = t12 * t48 + t45 * t50;
	t3 = t12 * t45 - t48 * t50;
	t2 = t10 * t48 + t45 * t51;
	t1 = t10 * t45 - t48 * t51;
	t5 = [0, t18 * pkin(3) - t57 * pkin(2) + t31 * t79 + t78 * (t18 * t45 - t31 * t75) + t59 * (t18 * t48 + t31 * t76) + t58 * (t31 * t61 - t57 * t46), t54 * t11 + t58 * t12, -t59 * t3 + t78 * t4, t3, (t11 * t38 - t4 * t37) * r_i_i_C(1) + (-t11 * t37 - t4 * t38) * r_i_i_C(2); 0, t16 * pkin(3) - t56 * pkin(2) + t30 * t79 + t78 * (t16 * t45 - t30 * t75) + t59 * (t16 * t48 + t30 * t76) + t58 * (t30 * t61 - t56 * t46), t58 * t10 + t54 * t9, -t59 * t1 + t78 * t2, t1, (-t2 * t37 + t9 * t38) * r_i_i_C(1) + (-t2 * t38 - t9 * t37) * r_i_i_C(2); 1, t28 * pkin(3) + t59 * (t28 * t48 + t45 * t68) + t78 * (t28 * t45 - t48 * t68) + (t49 * pkin(2) + t47 * t79 + t58 * (t46 * t49 + t47 * t61)) * t42, t54 * t23 + t58 * t24, -t59 * t13 + t78 * t14, t13, (-t14 * t37 + t23 * t38) * r_i_i_C(1) + (-t14 * t38 - t23 * t37) * r_i_i_C(2);];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end