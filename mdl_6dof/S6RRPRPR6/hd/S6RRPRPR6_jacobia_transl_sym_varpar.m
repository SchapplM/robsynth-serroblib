% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
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
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
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
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
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
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
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
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (204->42), mult. (503->67), div. (0->0), fcn. (658->10), ass. (0->35)
	t45 = pkin(4) - r_i_i_C(2);
	t44 = pkin(9) + r_i_i_C(1);
	t33 = cos(qJ(2));
	t43 = t33 * pkin(2);
	t28 = cos(pkin(6));
	t42 = t28 * t33;
	t26 = sin(pkin(6));
	t31 = sin(qJ(1));
	t41 = t31 * t26;
	t34 = cos(qJ(1));
	t40 = t34 * t26;
	t39 = r_i_i_C(3) + qJ(5);
	t25 = sin(pkin(11));
	t27 = cos(pkin(11));
	t30 = sin(qJ(2));
	t37 = t33 * t25 + t30 * t27;
	t18 = t37 * t28;
	t20 = t30 * t25 - t33 * t27;
	t9 = t34 * t18 - t31 * t20;
	t38 = t31 * t18 + t34 * t20;
	t29 = sin(qJ(4));
	t32 = cos(qJ(4));
	t1 = t9 * t29 + t32 * t40;
	t36 = t29 * t40 - t9 * t32;
	t35 = t39 * t29 + t45 * t32 + pkin(3);
	t24 = pkin(1) + t43;
	t19 = t28 * t30 * pkin(2) + (-pkin(8) - qJ(3)) * t26;
	t17 = t20 * t28;
	t16 = t37 * t26;
	t13 = t16 * t29 - t28 * t32;
	t11 = t31 * t17 - t34 * t37;
	t8 = -t34 * t17 - t31 * t37;
	t6 = t29 * t41 - t32 * t38;
	t5 = -t29 * t38 - t32 * t41;
	t2 = [-t9 * pkin(3) - t39 * t1 - t34 * t19 - t31 * t24 + t45 * t36 + t44 * t8, -t44 * t38 + (-t30 * t34 - t31 * t42) * pkin(2) + t35 * t11, t41, t39 * t6 - t45 * t5, t5, 0; -pkin(3) * t38 - t44 * t11 - t31 * t19 + t34 * t24 + t39 * t5 + t45 * t6, t44 * t9 + (-t30 * t31 + t34 * t42) * pkin(2) + t35 * t8, -t40, -t45 * t1 - t39 * t36, t1, 0; 0, t44 * t16 + (-t20 * t35 + t43) * t26, t28, t39 * (t16 * t32 + t28 * t29) - t45 * t13, t13, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (349->56), mult. (877->92), div. (0->0), fcn. (1158->12), ass. (0->41)
	t30 = sin(pkin(11));
	t35 = sin(qJ(2));
	t39 = cos(qJ(2));
	t50 = cos(pkin(11));
	t43 = -t35 * t30 + t39 * t50;
	t34 = sin(qJ(4));
	t38 = cos(qJ(4));
	t33 = sin(qJ(6));
	t37 = cos(qJ(6));
	t46 = t33 * r_i_i_C(1) + t37 * r_i_i_C(2) + qJ(5);
	t49 = r_i_i_C(3) + pkin(10) + pkin(4);
	t41 = t46 * t34 + t49 * t38 + pkin(3);
	t56 = -pkin(5) - pkin(9);
	t55 = t39 * pkin(2);
	t32 = cos(pkin(6));
	t54 = t32 * t39;
	t31 = sin(pkin(6));
	t36 = sin(qJ(1));
	t52 = t36 * t31;
	t40 = cos(qJ(1));
	t51 = t40 * t31;
	t24 = -t39 * t30 - t35 * t50;
	t21 = t24 * t32;
	t11 = -t40 * t21 + t36 * t43;
	t47 = -t36 * t21 - t40 * t43;
	t3 = t11 * t34 + t38 * t51;
	t45 = -t11 * t38 + t34 * t51;
	t44 = t37 * r_i_i_C(1) - t33 * r_i_i_C(2) - t56;
	t42 = t43 * t32;
	t29 = pkin(1) + t55;
	t22 = t32 * t35 * pkin(2) + (-pkin(8) - qJ(3)) * t31;
	t20 = t24 * t31;
	t19 = t43 * t31;
	t15 = -t20 * t34 - t32 * t38;
	t13 = t40 * t24 - t36 * t42;
	t10 = t36 * t24 + t40 * t42;
	t8 = t34 * t52 - t38 * t47;
	t7 = -t34 * t47 - t38 * t52;
	t2 = -t13 * t37 + t7 * t33;
	t1 = t13 * t33 + t7 * t37;
	t4 = [-t11 * pkin(3) + t44 * t10 - t40 * t22 - t36 * t29 - t46 * t3 + t49 * t45, (-t40 * t35 - t36 * t54) * pkin(2) - t44 * t47 + t41 * t13, t52, t46 * t8 - t49 * t7, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); -pkin(3) * t47 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(5) + t56 * t13 - t36 * t22 + t40 * t29 + t49 * t8, (-t36 * t35 + t40 * t54) * pkin(2) + t44 * t11 + t41 * t10, -t51, -t49 * t3 - t46 * t45, t3, (t10 * t33 + t3 * t37) * r_i_i_C(1) + (t10 * t37 - t3 * t33) * r_i_i_C(2); 0, t41 * t19 - t44 * t20 + t31 * t55, t32, t46 * (-t20 * t38 + t32 * t34) - t49 * t15, t15, (t15 * t37 + t19 * t33) * r_i_i_C(1) + (-t15 * t33 + t19 * t37) * r_i_i_C(2);];
	Ja_transl = t4;
else
	Ja_transl=NaN(3,6);
end