% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:55
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
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
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
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
	t11 = (pkin(8) + r_i_i_C(3)) * t5;
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
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (49->24), mult. (114->41), div. (0->0), fcn. (143->8), ass. (0->19)
	t16 = cos(qJ(2));
	t10 = sin(pkin(12));
	t12 = cos(pkin(12));
	t14 = sin(qJ(2));
	t6 = t10 * t14 - t16 * t12;
	t24 = -t16 * pkin(2) + t6 * r_i_i_C(1);
	t13 = cos(pkin(6));
	t21 = t13 * t16;
	t20 = -pkin(1) + t24;
	t19 = t10 * t16 + t12 * t14;
	t11 = sin(pkin(6));
	t4 = t19 * t13;
	t18 = -pkin(2) * t13 * t14 - t4 * r_i_i_C(1) + (r_i_i_C(3) + pkin(8) + qJ(3)) * t11;
	t17 = cos(qJ(1));
	t15 = sin(qJ(1));
	t3 = t6 * t13;
	t2 = t15 * t3 - t17 * t19;
	t1 = -t15 * t19 - t17 * t3;
	t5 = [-t1 * r_i_i_C(2) + t20 * t15 + t18 * t17, t2 * r_i_i_C(1) + (t15 * t4 + t17 * t6) * r_i_i_C(2) + (-t14 * t17 - t15 * t21) * pkin(2), t15 * t11, 0, 0, 0; t2 * r_i_i_C(2) + t18 * t15 - t20 * t17, t1 * r_i_i_C(1) + (t15 * t6 - t17 * t4) * r_i_i_C(2) + (-t14 * t15 + t17 * t21) * pkin(2), -t17 * t11, 0, 0, 0; 0, (-t19 * r_i_i_C(2) - t24) * t11, t13, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (128->40), mult. (313->67), div. (0->0), fcn. (405->10), ass. (0->32)
	t38 = r_i_i_C(3) + pkin(9);
	t27 = cos(qJ(2));
	t37 = t27 * pkin(2);
	t22 = cos(pkin(6));
	t36 = t22 * t27;
	t20 = sin(pkin(6));
	t25 = sin(qJ(1));
	t35 = t25 * t20;
	t28 = cos(qJ(1));
	t34 = t28 * t20;
	t23 = sin(qJ(4));
	t26 = cos(qJ(4));
	t19 = sin(pkin(12));
	t21 = cos(pkin(12));
	t24 = sin(qJ(2));
	t31 = t27 * t19 + t24 * t21;
	t12 = t31 * t22;
	t14 = t24 * t19 - t27 * t21;
	t5 = t28 * t12 - t25 * t14;
	t33 = t23 * t34 - t5 * t26;
	t32 = t25 * t12 + t28 * t14;
	t30 = t26 * r_i_i_C(1) - t23 * r_i_i_C(2) + pkin(3);
	t29 = t5 * t23 + t26 * t34;
	t18 = pkin(1) + t37;
	t13 = t22 * t24 * pkin(2) + (-pkin(8) - qJ(3)) * t20;
	t11 = t14 * t22;
	t10 = t31 * t20;
	t7 = t25 * t11 - t28 * t31;
	t4 = -t28 * t11 - t25 * t31;
	t2 = t23 * t35 - t26 * t32;
	t1 = t23 * t32 + t26 * t35;
	t3 = [-t5 * pkin(3) + t33 * r_i_i_C(1) + t29 * r_i_i_C(2) - t28 * t13 - t25 * t18 + t38 * t4, -t38 * t32 + (-t24 * t28 - t25 * t36) * pkin(2) + t30 * t7, t35, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; -pkin(3) * t32 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t25 * t13 + t28 * t18 - t38 * t7, t38 * t5 + (-t24 * t25 + t28 * t36) * pkin(2) + t30 * t4, -t34, -t29 * r_i_i_C(1) + t33 * r_i_i_C(2), 0, 0; 0, t38 * t10 + (-t14 * t30 + t37) * t20, t22, (-t10 * t23 + t22 * t26) * r_i_i_C(1) + (-t10 * t26 - t22 * t23) * r_i_i_C(2), 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (214->50), mult. (420->80), div. (0->0), fcn. (539->12), ass. (0->39)
	t28 = qJ(4) + qJ(5);
	t26 = sin(t28);
	t30 = sin(pkin(6));
	t38 = cos(qJ(1));
	t44 = t38 * t30;
	t22 = t26 * t44;
	t27 = cos(t28);
	t32 = cos(pkin(6));
	t29 = sin(pkin(12));
	t31 = cos(pkin(12));
	t34 = sin(qJ(2));
	t37 = cos(qJ(2));
	t41 = t37 * t29 + t34 * t31;
	t18 = t41 * t32;
	t20 = t34 * t29 - t37 * t31;
	t35 = sin(qJ(1));
	t9 = t38 * t18 - t35 * t20;
	t51 = (-t9 * t26 - t27 * t44) * r_i_i_C(1) + (-t9 * t27 + t22) * r_i_i_C(2);
	t42 = t35 * t18 + t38 * t20;
	t45 = t35 * t30;
	t5 = t26 * t42 + t27 * t45;
	t6 = t26 * t45 - t27 * t42;
	t50 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t33 = sin(qJ(4));
	t49 = pkin(4) * t33;
	t48 = t37 * pkin(2);
	t47 = r_i_i_C(3) + pkin(10) + pkin(9);
	t46 = t32 * t37;
	t16 = t41 * t30;
	t43 = (-t16 * t26 + t32 * t27) * r_i_i_C(1) + (-t16 * t27 - t32 * t26) * r_i_i_C(2);
	t36 = cos(qJ(4));
	t24 = t36 * pkin(4) + pkin(3);
	t40 = t27 * r_i_i_C(1) - t26 * r_i_i_C(2) + t24;
	t25 = pkin(1) + t48;
	t19 = t32 * t34 * pkin(2) + (-pkin(8) - qJ(3)) * t30;
	t17 = t20 * t32;
	t11 = t35 * t17 - t38 * t41;
	t8 = -t38 * t17 - t35 * t41;
	t1 = [t22 * r_i_i_C(1) - t35 * t25 - t40 * t9 + t47 * t8 + (-t19 + (r_i_i_C(2) * t27 + t49) * t30) * t38, -t47 * t42 + (-t34 * t38 - t35 * t46) * pkin(2) + t40 * t11, t45, (t33 * t42 + t36 * t45) * pkin(4) + t50, t50, 0; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t42 * t24 + t38 * t25 + (t30 * t49 - t19) * t35 - t47 * t11, t47 * t9 + (-t34 * t35 + t38 * t46) * pkin(2) + t40 * t8, -t44, (-t33 * t9 - t36 * t44) * pkin(4) + t51, t51, 0; 0, t47 * t16 + (-t20 * t40 + t48) * t30, t32, (-t16 * t33 + t32 * t36) * pkin(4) + t43, t43, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (482->67), mult. (933->105), div. (0->0), fcn. (1220->14), ass. (0->50)
	t76 = pkin(11) + r_i_i_C(3);
	t46 = sin(qJ(6));
	t50 = cos(qJ(6));
	t79 = t50 * r_i_i_C(1) - t46 * r_i_i_C(2) + pkin(5);
	t43 = sin(pkin(12));
	t48 = sin(qJ(2));
	t52 = cos(qJ(2));
	t66 = cos(pkin(12));
	t33 = -t52 * t43 - t48 * t66;
	t45 = cos(pkin(6));
	t30 = t33 * t45;
	t49 = sin(qJ(1));
	t53 = cos(qJ(1));
	t60 = -t48 * t43 + t52 * t66;
	t17 = -t53 * t30 + t49 * t60;
	t42 = qJ(4) + qJ(5);
	t40 = sin(t42);
	t41 = cos(t42);
	t44 = sin(pkin(6));
	t67 = t53 * t44;
	t10 = t17 * t41 - t40 * t67;
	t57 = t60 * t45;
	t16 = t49 * t33 + t53 * t57;
	t78 = t10 * t46 + t16 * t50;
	t77 = -t10 * t50 + t16 * t46;
	t51 = cos(qJ(4));
	t38 = t51 * pkin(4) + pkin(3);
	t55 = t76 * t40 + t79 * t41 + t38;
	t73 = t52 * pkin(2);
	t70 = t45 * t52;
	t68 = t49 * t44;
	t47 = sin(qJ(4));
	t64 = -t45 * t48 * pkin(2) + (t47 * pkin(4) + pkin(8) + qJ(3)) * t44;
	t63 = -t49 * t30 - t53 * t60;
	t54 = -pkin(10) - pkin(9);
	t61 = t46 * r_i_i_C(1) + t50 * r_i_i_C(2) - t54;
	t9 = -t17 * t40 - t41 * t67;
	t59 = t76 * t10 + t79 * t9;
	t13 = -t40 * t63 - t41 * t68;
	t14 = t40 * t68 - t41 * t63;
	t58 = -t79 * t13 + t76 * t14;
	t29 = t33 * t44;
	t25 = -t29 * t41 + t45 * t40;
	t56 = t76 * t25 + t79 * (t29 * t40 + t45 * t41);
	t39 = pkin(1) + t73;
	t28 = t60 * t44;
	t19 = t53 * t33 - t49 * t57;
	t2 = t14 * t50 - t19 * t46;
	t1 = -t14 * t46 - t19 * t50;
	t3 = [-t10 * pkin(5) + t77 * r_i_i_C(1) + t78 * r_i_i_C(2) - t16 * t54 - t17 * t38 - t49 * t39 + t64 * t53 + t76 * t9, (-t53 * t48 - t49 * t70) * pkin(2) - t61 * t63 + t55 * t19, t68, (t47 * t63 + t51 * t68) * pkin(4) + t58, t58, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t14 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t76 * t13 + t19 * t54 - t38 * t63 + t53 * t39 + t64 * t49, (-t49 * t48 + t53 * t70) * pkin(2) + t61 * t17 + t55 * t16, -t67, (-t17 * t47 - t51 * t67) * pkin(4) + t59, t59, -t78 * r_i_i_C(1) + t77 * r_i_i_C(2); 0, t55 * t28 - t61 * t29 + t44 * t73, t45, (t29 * t47 + t45 * t51) * pkin(4) + t56, t56, (-t25 * t46 - t28 * t50) * r_i_i_C(1) + (-t25 * t50 + t28 * t46) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end