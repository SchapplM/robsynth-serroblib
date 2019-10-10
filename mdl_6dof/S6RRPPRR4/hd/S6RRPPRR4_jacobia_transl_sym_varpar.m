% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:41
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
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
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
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
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
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
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (92->29), mult. (216->43), div. (0->0), fcn. (281->8), ass. (0->24)
	t32 = pkin(3) - r_i_i_C(2);
	t22 = cos(qJ(2));
	t31 = pkin(2) * t22;
	t19 = cos(pkin(6));
	t30 = t19 * t22;
	t29 = r_i_i_C(3) + qJ(4);
	t17 = sin(pkin(6));
	t20 = sin(qJ(2));
	t28 = -pkin(2) * t19 * t20 + (r_i_i_C(1) + pkin(8) + qJ(3)) * t17;
	t16 = sin(pkin(11));
	t18 = cos(pkin(11));
	t25 = t16 * t22 + t18 * t20;
	t10 = t25 * t19;
	t12 = t16 * t20 - t22 * t18;
	t21 = sin(qJ(1));
	t23 = cos(qJ(1));
	t27 = t10 * t23 - t21 * t12;
	t26 = t21 * t10 + t12 * t23;
	t24 = t12 * t19;
	t15 = pkin(1) + t31;
	t8 = t12 * t17;
	t5 = t21 * t24 - t23 * t25;
	t2 = -t21 * t25 - t23 * t24;
	t1 = [-t21 * t15 + t29 * t2 + t28 * t23 - t32 * t27, t32 * t5 - t29 * t26 + (-t20 * t23 - t21 * t30) * pkin(2), t21 * t17, -t5, 0, 0; t15 * t23 + t28 * t21 - t32 * t26 - t29 * t5, t32 * t2 + t29 * t27 + (-t20 * t21 + t23 * t30) * pkin(2), -t23 * t17, -t2, 0, 0; 0, -t32 * t8 + (t25 * t29 + t31) * t17, t19, t8, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (155->44), mult. (378->69), div. (0->0), fcn. (493->10), ass. (0->31)
	t21 = sin(pkin(11));
	t23 = cos(pkin(11));
	t26 = sin(qJ(2));
	t29 = cos(qJ(2));
	t15 = t26 * t21 - t29 * t23;
	t42 = t29 * pkin(2);
	t24 = cos(pkin(6));
	t41 = t24 * t29;
	t22 = sin(pkin(6));
	t27 = sin(qJ(1));
	t39 = t27 * t22;
	t30 = cos(qJ(1));
	t37 = t30 * t22;
	t36 = -r_i_i_C(3) - pkin(9) - pkin(3);
	t33 = t29 * t21 + t26 * t23;
	t13 = t33 * t24;
	t35 = t30 * t13 - t27 * t15;
	t34 = t27 * t13 + t30 * t15;
	t25 = sin(qJ(5));
	t28 = cos(qJ(5));
	t32 = t25 * r_i_i_C(1) + t28 * r_i_i_C(2) + qJ(4);
	t31 = t15 * t24;
	t20 = pkin(1) + t42;
	t18 = t28 * t37;
	t14 = t24 * t26 * pkin(2) + (-pkin(8) - qJ(3)) * t22;
	t11 = t15 * t22;
	t7 = t27 * t31 - t30 * t33;
	t4 = -t27 * t33 - t30 * t31;
	t2 = -t7 * t25 + t28 * t39;
	t1 = -t25 * t39 - t7 * t28;
	t3 = [t18 * r_i_i_C(1) - t27 * t20 + t32 * t4 + (-t14 + (-t25 * r_i_i_C(2) + pkin(4)) * t22) * t30 + t36 * t35, (-t26 * t30 - t27 * t41) * pkin(2) - t36 * t7 - t32 * t34, t39, -t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t7 * qJ(4) + t30 * t20 + (t22 * pkin(4) - t14) * t27 + t36 * t34, (-t26 * t27 + t30 * t41) * pkin(2) - t36 * t4 + t32 * t35, -t37, -t4, (t25 * t37 - t4 * t28) * r_i_i_C(1) + (t4 * t25 + t18) * r_i_i_C(2), 0; 0, t36 * t11 + (t32 * t33 + t42) * t22, t24, t11, (t11 * t28 - t24 * t25) * r_i_i_C(1) + (-t11 * t25 - t24 * t28) * r_i_i_C(2), 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (315->60), mult. (792->92), div. (0->0), fcn. (1045->12), ass. (0->41)
	t31 = sin(pkin(11));
	t36 = sin(qJ(2));
	t40 = cos(qJ(2));
	t50 = cos(pkin(11));
	t25 = -t40 * t31 - t36 * t50;
	t43 = -t36 * t31 + t40 * t50;
	t37 = sin(qJ(1));
	t41 = cos(qJ(1));
	t33 = cos(pkin(6));
	t51 = t25 * t33;
	t16 = t37 * t51 + t41 * t43;
	t63 = -t37 * t43 + t41 * t51;
	t35 = sin(qJ(5));
	t39 = cos(qJ(5));
	t34 = sin(qJ(6));
	t38 = cos(qJ(6));
	t46 = t38 * r_i_i_C(1) - t34 * r_i_i_C(2) + pkin(5);
	t60 = pkin(10) + r_i_i_C(3);
	t62 = t46 * t35 - t60 * t39 + qJ(4);
	t61 = pkin(3) + pkin(9);
	t59 = t40 * pkin(2);
	t58 = t33 * t40;
	t32 = sin(pkin(6));
	t55 = t37 * t32;
	t52 = t41 * t32;
	t49 = -t33 * t36 * pkin(2) + (pkin(4) + pkin(8) + qJ(3)) * t32;
	t22 = t43 * t33;
	t12 = t41 * t22 + t37 * t25;
	t8 = t12 * t35 + t39 * t52;
	t45 = -t12 * t39 + t35 * t52;
	t44 = t34 * r_i_i_C(1) + t38 * r_i_i_C(2) + t61;
	t30 = pkin(1) + t59;
	t21 = t25 * t32;
	t20 = t43 * t32;
	t18 = -t20 * t35 + t33 * t39;
	t15 = -t37 * t22 + t41 * t25;
	t4 = -t15 * t35 + t39 * t55;
	t3 = t15 * t39 + t35 * t55;
	t2 = t16 * t34 + t4 * t38;
	t1 = t16 * t38 - t4 * t34;
	t5 = [t12 * qJ(4) - t37 * t30 + t49 * t41 + t44 * t63 + t60 * t45 + t46 * t8, (-t41 * t36 - t37 * t58) * pkin(2) + t44 * t15 + t62 * t16, t55, -t15, -t46 * t3 + t60 * t4, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t4 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t15 * qJ(4) + t61 * t16 + t60 * t3 + t41 * t30 + t49 * t37, (-t37 * t36 + t41 * t58) * pkin(2) + t44 * t12 - t62 * t63, -t52, -t12, t46 * t45 - t60 * t8, (t8 * t34 - t38 * t63) * r_i_i_C(1) + (t34 * t63 + t8 * t38) * r_i_i_C(2); 0, t44 * t20 - t62 * t21 + t32 * t59, t33, -t20, t60 * t18 + t46 * (-t20 * t39 - t33 * t35), (-t18 * t34 - t21 * t38) * r_i_i_C(1) + (-t18 * t38 + t21 * t34) * r_i_i_C(2);];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end