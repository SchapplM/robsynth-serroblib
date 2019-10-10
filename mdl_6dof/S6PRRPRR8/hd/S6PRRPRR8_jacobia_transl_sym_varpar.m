% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR8_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR8_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR8_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
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
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (111->35), mult. (307->66), div. (0->0), fcn. (394->10), ass. (0->33)
	t42 = pkin(3) - r_i_i_C(2);
	t19 = sin(pkin(7));
	t20 = sin(pkin(6));
	t41 = t19 * t20;
	t23 = cos(pkin(6));
	t40 = t19 * t23;
	t22 = cos(pkin(7));
	t24 = sin(qJ(3));
	t39 = t22 * t24;
	t26 = cos(qJ(3));
	t38 = t22 * t26;
	t25 = sin(qJ(2));
	t37 = t23 * t25;
	t27 = cos(qJ(2));
	t36 = t23 * t27;
	t35 = t24 * t25;
	t34 = t24 * t27;
	t33 = t25 * t26;
	t32 = t26 * t27;
	t31 = r_i_i_C(3) + qJ(4);
	t30 = (pkin(9) + r_i_i_C(1)) * t19;
	t18 = sin(pkin(12));
	t21 = cos(pkin(12));
	t13 = -t18 * t25 + t21 * t36;
	t29 = -t13 * t22 + t21 * t41;
	t15 = -t18 * t36 - t21 * t25;
	t28 = t15 * t22 + t18 * t41;
	t16 = -t18 * t37 + t21 * t27;
	t14 = t18 * t27 + t21 * t37;
	t9 = -t26 * t40 + (-t22 * t32 + t35) * t20;
	t3 = t16 * t24 - t28 * t26;
	t1 = t14 * t24 + t29 * t26;
	t2 = [0, t15 * pkin(2) + t42 * (t15 * t26 - t16 * t39) + t31 * (t15 * t24 + t16 * t38) + t16 * t30, t31 * (t16 * t26 + t28 * t24) - t42 * t3, t3, 0, 0; 0, t13 * pkin(2) + t42 * (t13 * t26 - t14 * t39) + t31 * (t13 * t24 + t14 * t38) + t14 * t30, t31 * (t14 * t26 - t29 * t24) - t42 * t1, t1, 0, 0; 1, (t42 * (-t22 * t35 + t32) + t31 * (t22 * t33 + t34) + pkin(2) * t27 + t25 * t30) * t20, -t42 * t9 + t31 * (t24 * t40 + (t22 * t34 + t33) * t20), t9, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (206->52), mult. (575->98), div. (0->0), fcn. (742->12), ass. (0->40)
	t26 = sin(pkin(7));
	t27 = sin(pkin(6));
	t54 = t26 * t27;
	t29 = cos(pkin(7));
	t53 = t27 * t29;
	t31 = sin(qJ(3));
	t52 = t29 * t31;
	t34 = cos(qJ(3));
	t51 = t29 * t34;
	t32 = sin(qJ(2));
	t50 = t31 * t32;
	t35 = cos(qJ(2));
	t49 = t31 * t35;
	t48 = t32 * t34;
	t47 = t34 * t35;
	t46 = cos(pkin(6));
	t45 = sin(pkin(12));
	t44 = r_i_i_C(3) + pkin(10) + pkin(3);
	t28 = cos(pkin(12));
	t43 = t28 * t54;
	t42 = t27 * t45;
	t41 = t28 * t46;
	t40 = t46 * t26;
	t39 = t26 * t42;
	t38 = t46 * t45;
	t30 = sin(qJ(5));
	t33 = cos(qJ(5));
	t37 = t30 * r_i_i_C(1) + t33 * r_i_i_C(2) + qJ(4);
	t36 = (t33 * r_i_i_C(1) - t30 * r_i_i_C(2) + pkin(4) + pkin(9)) * t26;
	t21 = t28 * t35 - t32 * t38;
	t20 = -t28 * t32 - t35 * t38;
	t19 = t32 * t41 + t45 * t35;
	t18 = -t45 * t32 + t35 * t41;
	t17 = t46 * t29 - t35 * t54;
	t12 = -t20 * t26 + t29 * t42;
	t11 = -t18 * t26 - t28 * t53;
	t9 = t27 * t50 - t34 * t40 - t47 * t53;
	t3 = -t20 * t51 + t21 * t31 - t34 * t39;
	t1 = -t18 * t51 + t19 * t31 + t34 * t43;
	t2 = [0, t20 * pkin(2) + t37 * (t20 * t31 + t21 * t51) + t21 * t36 + t44 * (t20 * t34 - t21 * t52), t37 * (t21 * t34 + (t20 * t29 + t39) * t31) - t44 * t3, t3, (-t12 * t30 + t3 * t33) * r_i_i_C(1) + (-t12 * t33 - t3 * t30) * r_i_i_C(2), 0; 0, t18 * pkin(2) + t37 * (t18 * t31 + t19 * t51) + t19 * t36 + t44 * (t18 * t34 - t19 * t52), t37 * (t19 * t34 + (t18 * t29 - t43) * t31) - t44 * t1, t1, (t1 * t33 - t11 * t30) * r_i_i_C(1) + (-t1 * t30 - t11 * t33) * r_i_i_C(2), 0; 1, (t37 * (t29 * t48 + t49) + t35 * pkin(2) + t32 * t36 + t44 * (-t29 * t50 + t47)) * t27, -t44 * t9 + t37 * (t31 * t40 + (t29 * t49 + t48) * t27), t9, (-t17 * t30 + t9 * t33) * r_i_i_C(1) + (-t17 * t33 - t9 * t30) * r_i_i_C(2), 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:10
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (441->76), mult. (1239->140), div. (0->0), fcn. (1617->14), ass. (0->56)
	t73 = r_i_i_C(3) + pkin(11);
	t38 = sin(pkin(7));
	t39 = sin(pkin(6));
	t72 = t38 * t39;
	t43 = sin(qJ(5));
	t71 = t38 * t43;
	t47 = cos(qJ(5));
	t70 = t38 * t47;
	t41 = cos(pkin(7));
	t69 = t39 * t41;
	t44 = sin(qJ(3));
	t68 = t41 * t44;
	t48 = cos(qJ(3));
	t67 = t41 * t48;
	t45 = sin(qJ(2));
	t66 = t44 * t45;
	t49 = cos(qJ(2));
	t65 = t44 * t49;
	t64 = t45 * t48;
	t63 = t48 * t49;
	t62 = cos(pkin(6));
	t61 = sin(pkin(12));
	t40 = cos(pkin(12));
	t60 = t40 * t72;
	t59 = t45 * t72;
	t58 = (pkin(4) + pkin(9)) * t38;
	t57 = t39 * t61;
	t56 = t40 * t62;
	t55 = t62 * t38;
	t54 = t38 * t57;
	t53 = t62 * t61;
	t42 = sin(qJ(6));
	t46 = cos(qJ(6));
	t52 = t46 * r_i_i_C(1) - t42 * r_i_i_C(2) + pkin(5);
	t51 = t42 * r_i_i_C(1) + t46 * r_i_i_C(2) + pkin(3) + pkin(10);
	t50 = t52 * t43 - t73 * t47 + qJ(4);
	t33 = t40 * t49 - t45 * t53;
	t32 = -t40 * t45 - t49 * t53;
	t31 = t45 * t56 + t61 * t49;
	t30 = -t61 * t45 + t49 * t56;
	t29 = t62 * t41 - t49 * t72;
	t27 = (t41 * t64 + t65) * t39;
	t24 = -t32 * t38 + t41 * t57;
	t23 = -t30 * t38 - t40 * t69;
	t22 = t44 * t55 + (t41 * t65 + t64) * t39;
	t21 = t39 * t66 - t48 * t55 - t63 * t69;
	t17 = t32 * t44 + t33 * t67;
	t15 = t30 * t44 + t31 * t67;
	t14 = t21 * t43 + t29 * t47;
	t12 = t33 * t48 + (t32 * t41 + t54) * t44;
	t11 = -t32 * t67 + t33 * t44 - t48 * t54;
	t10 = t31 * t48 + (t30 * t41 - t60) * t44;
	t9 = -t30 * t67 + t31 * t44 + t48 * t60;
	t4 = t11 * t43 + t24 * t47;
	t2 = t23 * t47 + t9 * t43;
	t1 = [0, t32 * pkin(2) + t17 * qJ(4) - t73 * (t17 * t47 - t33 * t71) + t33 * t58 + t52 * (t17 * t43 + t33 * t70) + t51 * (t32 * t48 - t33 * t68), -t51 * t11 + t50 * t12, t11, t73 * t4 + t52 * (t11 * t47 - t24 * t43), (t12 * t46 - t4 * t42) * r_i_i_C(1) + (-t12 * t42 - t4 * t46) * r_i_i_C(2); 0, t30 * pkin(2) + t15 * qJ(4) - t73 * (t15 * t47 - t31 * t71) + t31 * t58 + t52 * (t15 * t43 + t31 * t70) + t51 * (t30 * t48 - t31 * t68), t50 * t10 - t51 * t9, t9, t73 * t2 + t52 * (-t23 * t43 + t9 * t47), (t10 * t46 - t2 * t42) * r_i_i_C(1) + (-t10 * t42 - t2 * t46) * r_i_i_C(2); 1, t27 * qJ(4) + t52 * (t27 * t43 + t47 * t59) + t73 * (-t27 * t47 + t43 * t59) + (t51 * (-t41 * t66 + t63) + t49 * pkin(2) + t45 * t58) * t39, -t51 * t21 + t50 * t22, t21, t73 * t14 + t52 * (t21 * t47 - t29 * t43), (-t14 * t42 + t22 * t46) * r_i_i_C(1) + (-t14 * t46 - t22 * t42) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end