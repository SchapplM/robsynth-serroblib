% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:09
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
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
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
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
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (49->24), mult. (114->41), div. (0->0), fcn. (143->8), ass. (0->19)
	t16 = cos(qJ(2));
	t10 = sin(pkin(11));
	t12 = cos(pkin(11));
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
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
	% DurationCPUTime: 0.20s
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
	t19 = sin(pkin(11));
	t21 = cos(pkin(11));
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
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (186->52), mult. (378->81), div. (0->0), fcn. (487->12), ass. (0->37)
	t24 = sin(pkin(11));
	t26 = cos(pkin(11));
	t30 = sin(qJ(2));
	t33 = cos(qJ(2));
	t14 = t30 * t24 - t33 * t26;
	t29 = sin(qJ(4));
	t46 = t29 * pkin(4);
	t45 = t33 * pkin(2);
	t44 = r_i_i_C(3) + qJ(5) + pkin(9);
	t27 = cos(pkin(6));
	t43 = t27 * t33;
	t25 = sin(pkin(6));
	t31 = sin(qJ(1));
	t41 = t31 * t25;
	t34 = cos(qJ(1));
	t39 = t34 * t25;
	t37 = t33 * t24 + t30 * t26;
	t12 = t37 * t27;
	t5 = t34 * t12 - t31 * t14;
	t38 = t31 * t12 + t34 * t14;
	t32 = cos(qJ(4));
	t19 = t32 * pkin(4) + pkin(3);
	t23 = qJ(4) + pkin(12);
	t21 = sin(t23);
	t22 = cos(t23);
	t36 = t22 * r_i_i_C(1) - t21 * r_i_i_C(2) + t19;
	t35 = t14 * t27;
	t20 = pkin(1) + t45;
	t16 = t21 * t39;
	t13 = t27 * t30 * pkin(2) + (-pkin(8) - qJ(3)) * t25;
	t11 = t37 * t25;
	t10 = t14 * t25;
	t7 = t31 * t35 - t34 * t37;
	t4 = -t31 * t37 - t34 * t35;
	t2 = t21 * t41 - t22 * t38;
	t1 = t21 * t38 + t22 * t41;
	t3 = [t16 * r_i_i_C(1) - t31 * t20 - t36 * t5 + t44 * t4 + (-t13 + (t22 * r_i_i_C(2) + t46) * t25) * t34, -t44 * t38 + (-t30 * t34 - t31 * t43) * pkin(2) + t36 * t7, t41, t1 * r_i_i_C(1) - t2 * r_i_i_C(2) + (t29 * t38 + t32 * t41) * pkin(4), -t7, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t38 * t19 + t34 * t20 - t44 * t7 + (t25 * t46 - t13) * t31, t44 * t5 + (-t30 * t31 + t34 * t43) * pkin(2) + t36 * t4, -t39, (-t5 * t21 - t22 * t39) * r_i_i_C(1) + (-t5 * t22 + t16) * r_i_i_C(2) + (-t5 * t29 - t32 * t39) * pkin(4), -t4, 0; 0, -t36 * t10 + t44 * t11 + t25 * t45, t27, (-t11 * t21 + t27 * t22) * r_i_i_C(1) + (-t11 * t22 - t27 * t21) * r_i_i_C(2) + (-t11 * t29 + t27 * t32) * pkin(4), t10, 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:42
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (397->70), mult. (792->105), div. (0->0), fcn. (1039->14), ass. (0->47)
	t34 = sin(pkin(11));
	t40 = sin(qJ(2));
	t44 = cos(qJ(2));
	t55 = cos(pkin(11));
	t24 = -t44 * t34 - t40 * t55;
	t41 = sin(qJ(1));
	t45 = cos(qJ(1));
	t36 = cos(pkin(6));
	t48 = -t40 * t34 + t44 * t55;
	t47 = t48 * t36;
	t10 = t41 * t24 + t45 * t47;
	t38 = sin(qJ(6));
	t21 = t24 * t36;
	t11 = -t45 * t21 + t41 * t48;
	t33 = qJ(4) + pkin(12);
	t31 = sin(t33);
	t32 = cos(t33);
	t35 = sin(pkin(6));
	t56 = t45 * t35;
	t4 = t11 * t32 - t31 * t56;
	t42 = cos(qJ(6));
	t65 = t10 * t42 + t4 * t38;
	t64 = t10 * t38 - t4 * t42;
	t43 = cos(qJ(4));
	t29 = t43 * pkin(4) + pkin(3);
	t51 = t42 * r_i_i_C(1) - t38 * r_i_i_C(2) + pkin(5);
	t63 = pkin(10) + r_i_i_C(3);
	t46 = t63 * t31 + t51 * t32 + t29;
	t62 = t44 * pkin(2);
	t59 = t36 * t44;
	t57 = t41 * t35;
	t39 = sin(qJ(4));
	t53 = -t36 * t40 * pkin(2) + (t39 * pkin(4) + pkin(8) + qJ(3)) * t35;
	t52 = -t41 * t21 - t45 * t48;
	t37 = -qJ(5) - pkin(9);
	t50 = t38 * r_i_i_C(1) + t42 * r_i_i_C(2) - t37;
	t49 = -t11 * t31 - t32 * t56;
	t30 = pkin(1) + t62;
	t20 = t24 * t35;
	t19 = t48 * t35;
	t16 = -t20 * t32 + t36 * t31;
	t13 = t45 * t24 - t41 * t47;
	t8 = t31 * t57 - t32 * t52;
	t7 = -t31 * t52 - t32 * t57;
	t2 = -t13 * t38 + t8 * t42;
	t1 = -t13 * t42 - t8 * t38;
	t3 = [-t4 * pkin(5) + t64 * r_i_i_C(1) + t65 * r_i_i_C(2) - t10 * t37 - t11 * t29 - t41 * t30 + t53 * t45 + t63 * t49, (-t45 * t40 - t41 * t59) * pkin(2) - t50 * t52 + t46 * t13, t57, t63 * t8 + (t39 * t52 + t43 * t57) * pkin(4) - t51 * t7, -t13, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t8 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t13 * t37 - t29 * t52 + t45 * t30 + t53 * t41 + t63 * t7, (-t41 * t40 + t45 * t59) * pkin(2) + t50 * t11 + t46 * t10, -t56, t63 * t4 + (-t11 * t39 - t43 * t56) * pkin(4) + t51 * t49, -t10, -t65 * r_i_i_C(1) + t64 * r_i_i_C(2); 0, t46 * t19 - t50 * t20 + t35 * t62, t36, t63 * t16 + (t20 * t39 + t36 * t43) * pkin(4) + t51 * (t20 * t31 + t36 * t32), -t19, (-t16 * t38 - t19 * t42) * r_i_i_C(1) + (-t16 * t42 + t19 * t38) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end