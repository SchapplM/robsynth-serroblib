% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:44
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP10_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP10_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP10_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:44:47
	% EndTime: 2019-10-10 10:44:48
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:44:47
	% EndTime: 2019-10-10 10:44:47
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
	% StartTime: 2019-10-10 10:44:47
	% EndTime: 2019-10-10 10:44:48
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
	% StartTime: 2019-10-10 10:44:47
	% EndTime: 2019-10-10 10:44:48
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (51->25), mult. (121->42), div. (0->0), fcn. (150->8), ass. (0->19)
	t10 = sin(qJ(1));
	t7 = sin(pkin(6));
	t19 = t10 * t7;
	t12 = cos(qJ(1));
	t18 = t12 * t7;
	t17 = r_i_i_C(3) + qJ(3);
	t16 = cos(pkin(6));
	t15 = t10 * t16;
	t14 = t12 * t16;
	t6 = sin(pkin(11));
	t8 = cos(pkin(11));
	t13 = t8 * r_i_i_C(1) - t6 * r_i_i_C(2) + pkin(2);
	t11 = cos(qJ(2));
	t9 = sin(qJ(2));
	t4 = t12 * t11 - t9 * t15;
	t3 = t11 * t15 + t12 * t9;
	t2 = t10 * t11 + t9 * t14;
	t1 = t10 * t9 - t11 * t14;
	t5 = [(t6 * t18 - t2 * t8) * r_i_i_C(1) + (t8 * t18 + t2 * t6) * r_i_i_C(2) - t2 * pkin(2) - t10 * pkin(1) + pkin(8) * t18 - t17 * t1, -t13 * t3 + t17 * t4, t3, 0, 0, 0; (t6 * t19 + t4 * t8) * r_i_i_C(1) + (t8 * t19 - t4 * t6) * r_i_i_C(2) + t4 * pkin(2) + t12 * pkin(1) + pkin(8) * t19 + t17 * t3, -t13 * t1 + t17 * t2, t1, 0, 0, 0; 0, (t13 * t11 + t17 * t9) * t7, -t7 * t11, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:44:47
	% EndTime: 2019-10-10 10:44:48
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (102->33), mult. (168->53), div. (0->0), fcn. (207->10), ass. (0->28)
	t30 = r_i_i_C(3) + pkin(9) + qJ(3);
	t14 = sin(pkin(6));
	t17 = sin(qJ(2));
	t29 = t14 * t17;
	t18 = sin(qJ(1));
	t28 = t14 * t18;
	t20 = cos(qJ(1));
	t27 = t14 * t20;
	t26 = t18 * t17;
	t19 = cos(qJ(2));
	t25 = t18 * t19;
	t24 = t20 * t17;
	t23 = t20 * t19;
	t22 = pkin(3) * sin(pkin(11)) + pkin(8);
	t12 = pkin(11) + qJ(4);
	t10 = sin(t12);
	t11 = cos(t12);
	t9 = cos(pkin(11)) * pkin(3) + pkin(2);
	t21 = t11 * r_i_i_C(1) - t10 * r_i_i_C(2) + t9;
	t15 = cos(pkin(6));
	t7 = t10 * t27;
	t6 = -t15 * t26 + t23;
	t5 = t15 * t25 + t24;
	t4 = t15 * t24 + t25;
	t3 = -t15 * t23 + t26;
	t2 = t10 * t28 + t6 * t11;
	t1 = -t6 * t10 + t11 * t28;
	t8 = [-t18 * pkin(1) + t7 * r_i_i_C(1) - t21 * t4 - t30 * t3 + (t11 * r_i_i_C(2) + t22) * t27, -t21 * t5 + t30 * t6, t5, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t20 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t22 * t28 + t30 * t5 + t6 * t9, -t21 * t3 + t30 * t4, t3, (-t4 * t10 - t11 * t27) * r_i_i_C(1) + (-t4 * t11 + t7) * r_i_i_C(2), 0, 0; 0, (t30 * t17 + t21 * t19) * t14, -t14 * t19, (-t10 * t29 + t15 * t11) * r_i_i_C(1) + (-t15 * t10 - t11 * t29) * r_i_i_C(2), 0, 0;];
	Ja_transl = t8;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:44:47
	% EndTime: 2019-10-10 10:44:48
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (242->51), mult. (399->82), div. (0->0), fcn. (505->12), ass. (0->38)
	t26 = sin(qJ(2));
	t27 = sin(qJ(1));
	t29 = cos(qJ(2));
	t30 = cos(qJ(1));
	t38 = cos(pkin(6));
	t36 = t30 * t38;
	t11 = t27 * t26 - t29 * t36;
	t25 = sin(qJ(5));
	t28 = cos(qJ(5));
	t12 = t26 * t36 + t27 * t29;
	t21 = pkin(11) + qJ(4);
	t19 = sin(t21);
	t20 = cos(t21);
	t23 = sin(pkin(6));
	t39 = t23 * t30;
	t4 = t12 * t20 - t19 * t39;
	t48 = -t11 * t28 + t4 * t25;
	t47 = -t11 * t25 - t4 * t28;
	t18 = cos(pkin(11)) * pkin(3) + pkin(2);
	t34 = t28 * r_i_i_C(1) - t25 * r_i_i_C(2) + pkin(4);
	t45 = pkin(10) + r_i_i_C(3);
	t46 = t45 * t19 + t34 * t20 + t18;
	t42 = t23 * t26;
	t41 = t23 * t27;
	t40 = t23 * t29;
	t37 = t27 * t38;
	t35 = t23 * (pkin(3) * sin(pkin(11)) + pkin(8));
	t24 = -pkin(9) - qJ(3);
	t33 = t25 * r_i_i_C(1) + t28 * r_i_i_C(2) - t24;
	t32 = -t12 * t19 - t20 * t39;
	t14 = -t26 * t37 + t30 * t29;
	t13 = t30 * t26 + t29 * t37;
	t10 = t38 * t19 + t20 * t42;
	t8 = t14 * t20 + t19 * t41;
	t7 = t14 * t19 - t20 * t41;
	t2 = t13 * t25 + t8 * t28;
	t1 = t13 * t28 - t8 * t25;
	t3 = [-t27 * pkin(1) - t4 * pkin(4) + t47 * r_i_i_C(1) + t48 * r_i_i_C(2) + t11 * t24 - t12 * t18 + t30 * t35 + t45 * t32, -t13 * t46 + t33 * t14, t13, -t34 * t7 + t45 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t30 * pkin(1) + t8 * pkin(4) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t13 * t24 + t14 * t18 + t27 * t35 + t45 * t7, -t11 * t46 + t33 * t12, t11, t34 * t32 + t45 * t4, -t48 * r_i_i_C(1) + t47 * r_i_i_C(2), 0; 0, (t33 * t26 + t29 * t46) * t23, -t40, t45 * t10 + t34 * (-t19 * t42 + t38 * t20), (-t10 * t25 - t28 * t40) * r_i_i_C(1) + (-t10 * t28 + t25 * t40) * r_i_i_C(2), 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:44:47
	% EndTime: 2019-10-10 10:44:48
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (378->63), mult. (631->102), div. (0->0), fcn. (810->12), ass. (0->43)
	t40 = sin(qJ(2));
	t41 = sin(qJ(1));
	t43 = cos(qJ(2));
	t44 = cos(qJ(1));
	t52 = cos(pkin(6));
	t49 = t44 * t52;
	t26 = t40 * t49 + t41 * t43;
	t35 = pkin(11) + qJ(4);
	t33 = sin(t35);
	t34 = cos(t35);
	t37 = sin(pkin(6));
	t56 = t37 * t44;
	t14 = t26 * t34 - t33 * t56;
	t25 = t41 * t40 - t43 * t49;
	t39 = sin(qJ(5));
	t42 = cos(qJ(5));
	t1 = t14 * t39 - t25 * t42;
	t66 = t14 * t42 + t25 * t39;
	t32 = cos(pkin(11)) * pkin(3) + pkin(2);
	t63 = pkin(10) + r_i_i_C(2);
	t65 = pkin(4) * t34 + t63 * t33 + t32;
	t53 = r_i_i_C(3) + qJ(6);
	t64 = pkin(5) + r_i_i_C(1);
	t45 = t53 * t39 + t64 * t42 + pkin(4);
	t60 = t34 * t39;
	t59 = t34 * t42;
	t58 = t37 * t40;
	t57 = t37 * t41;
	t55 = t39 * t43;
	t54 = t42 * t43;
	t50 = t41 * t52;
	t48 = t37 * (pkin(3) * sin(pkin(11)) + pkin(8));
	t47 = -t26 * t33 - t34 * t56;
	t38 = -pkin(9) - qJ(3);
	t28 = -t40 * t50 + t44 * t43;
	t27 = t44 * t40 + t43 * t50;
	t22 = t52 * t33 + t34 * t58;
	t18 = t28 * t34 + t33 * t57;
	t17 = t28 * t33 - t34 * t57;
	t11 = t22 * t39 + t37 * t54;
	t6 = t18 * t42 + t27 * t39;
	t5 = t18 * t39 - t27 * t42;
	t2 = [-t41 * pkin(1) - t14 * pkin(4) - t53 * t1 + t25 * t38 - t26 * t32 + t44 * t48 + t63 * t47 - t64 * t66, -t28 * t38 + t53 * (-t27 * t60 - t28 * t42) + t64 * (-t27 * t59 + t28 * t39) - t65 * t27, t27, -t45 * t17 + t63 * t18, -t64 * t5 + t53 * t6, t5; t44 * pkin(1) + t18 * pkin(4) + t63 * t17 - t27 * t38 + t28 * t32 + t41 * t48 + t53 * t5 + t64 * t6, -t26 * t38 + t64 * (-t25 * t59 + t26 * t39) + t53 * (-t25 * t60 - t26 * t42) - t65 * t25, t25, t63 * t14 + t45 * t47, -t64 * t1 + t53 * t66, t1; 0, (t64 * (t34 * t54 + t39 * t40) + t53 * (t34 * t55 - t40 * t42) - t38 * t40 + t65 * t43) * t37, -t37 * t43, t63 * t22 + t45 * (-t33 * t58 + t34 * t52), t53 * (t22 * t42 - t37 * t55) - t64 * t11, t11;];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,6);
end