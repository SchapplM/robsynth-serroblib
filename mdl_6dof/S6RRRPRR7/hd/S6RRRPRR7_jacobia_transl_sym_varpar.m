% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:04
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR7_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR7_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR7_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.11s
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
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (61->29), mult. (150->51), div. (0->0), fcn. (184->8), ass. (0->26)
	t27 = -r_i_i_C(3) - pkin(9);
	t13 = sin(qJ(1));
	t9 = sin(pkin(6));
	t26 = t13 * t9;
	t14 = cos(qJ(3));
	t25 = t14 * t9;
	t16 = cos(qJ(1));
	t24 = t16 * t9;
	t12 = sin(qJ(2));
	t23 = t12 * t13;
	t22 = t12 * t16;
	t15 = cos(qJ(2));
	t21 = t13 * t15;
	t20 = t15 * t16;
	t11 = sin(qJ(3));
	t10 = cos(pkin(6));
	t4 = t10 * t22 + t21;
	t19 = t11 * t24 - t14 * t4;
	t18 = t14 * r_i_i_C(1) - t11 * r_i_i_C(2) + pkin(2);
	t17 = t4 * t11 + t14 * t24;
	t6 = -t10 * t23 + t20;
	t5 = t10 * t21 + t22;
	t3 = -t10 * t20 + t23;
	t2 = t11 * t26 + t14 * t6;
	t1 = -t11 * t6 + t13 * t25;
	t7 = [-t13 * pkin(1) - t4 * pkin(2) + pkin(8) * t24 + t19 * r_i_i_C(1) + t17 * r_i_i_C(2) + t27 * t3, -t18 * t5 - t27 * t6, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; pkin(1) * t16 + t6 * pkin(2) + pkin(8) * t26 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t27 * t5, -t18 * t3 - t27 * t4, -t17 * r_i_i_C(1) + t19 * r_i_i_C(2), 0, 0, 0; 0, (-t27 * t12 + t18 * t15) * t9, (-t11 * t12 * t9 + t10 * t14) * r_i_i_C(1) + (-t10 * t11 - t12 * t25) * r_i_i_C(2), 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (110->39), mult. (189->62), div. (0->0), fcn. (230->10), ass. (0->30)
	t31 = r_i_i_C(3) + qJ(4) + pkin(9);
	t13 = sin(pkin(6));
	t17 = sin(qJ(2));
	t30 = t13 * t17;
	t18 = sin(qJ(1));
	t29 = t13 * t18;
	t21 = cos(qJ(1));
	t28 = t13 * t21;
	t27 = t17 * t18;
	t26 = t17 * t21;
	t20 = cos(qJ(2));
	t25 = t18 * t20;
	t24 = t20 * t21;
	t16 = sin(qJ(3));
	t23 = t16 * pkin(3) + pkin(8);
	t12 = qJ(3) + pkin(12);
	t10 = sin(t12);
	t11 = cos(t12);
	t19 = cos(qJ(3));
	t9 = t19 * pkin(3) + pkin(2);
	t22 = t11 * r_i_i_C(1) - t10 * r_i_i_C(2) + t9;
	t14 = cos(pkin(6));
	t7 = t10 * t28;
	t6 = -t14 * t27 + t24;
	t5 = t14 * t25 + t26;
	t4 = t14 * t26 + t25;
	t3 = -t14 * t24 + t27;
	t2 = t10 * t29 + t11 * t6;
	t1 = -t10 * t6 + t11 * t29;
	t8 = [-t18 * pkin(1) + t7 * r_i_i_C(1) - t22 * t4 - t31 * t3 + (t11 * r_i_i_C(2) + t23) * t28, -t22 * t5 + t31 * t6, r_i_i_C(1) * t1 - r_i_i_C(2) * t2 + (-t16 * t6 + t19 * t29) * pkin(3), t5, 0, 0; pkin(1) * t21 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t23 * t29 + t31 * t5 + t6 * t9, -t22 * t3 + t31 * t4, (-t4 * t10 - t11 * t28) * r_i_i_C(1) + (-t11 * t4 + t7) * r_i_i_C(2) + (-t4 * t16 - t19 * t28) * pkin(3), t3, 0, 0; 0, (t31 * t17 + t22 * t20) * t13, (-t10 * t30 + t11 * t14) * r_i_i_C(1) + (-t10 * t14 - t11 * t30) * r_i_i_C(2) + (t14 * t19 - t16 * t30) * pkin(3), -t13 * t20, 0, 0;];
	Ja_transl = t8;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (204->42), mult. (238->61), div. (0->0), fcn. (289->12), ass. (0->34)
	t26 = cos(pkin(6));
	t27 = sin(qJ(2));
	t30 = cos(qJ(1));
	t33 = t30 * t27;
	t28 = sin(qJ(1));
	t29 = cos(qJ(2));
	t34 = t28 * t29;
	t10 = t26 * t33 + t34;
	t24 = qJ(3) + pkin(12);
	t21 = qJ(5) + t24;
	t19 = sin(t21);
	t25 = sin(pkin(6));
	t36 = t25 * t30;
	t13 = t19 * t36;
	t20 = cos(t21);
	t43 = (-t10 * t19 - t20 * t36) * r_i_i_C(1) + (-t10 * t20 + t13) * r_i_i_C(2);
	t32 = t30 * t29;
	t35 = t28 * t27;
	t12 = -t26 * t35 + t32;
	t37 = t25 * t28;
	t5 = -t12 * t19 + t20 * t37;
	t6 = t12 * t20 + t19 * t37;
	t42 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t38 = t25 * t27;
	t41 = (-t19 * t38 + t26 * t20) * r_i_i_C(1) + (-t26 * t19 - t20 * t38) * r_i_i_C(2);
	t15 = pkin(4) * sin(t24) + sin(qJ(3)) * pkin(3);
	t40 = pkin(8) + t15;
	t39 = r_i_i_C(3) + pkin(10) + qJ(4) + pkin(9);
	t16 = pkin(4) * cos(t24) + cos(qJ(3)) * pkin(3);
	t14 = pkin(2) + t16;
	t31 = t20 * r_i_i_C(1) - t19 * r_i_i_C(2) + t14;
	t11 = t26 * t34 + t33;
	t9 = -t26 * t32 + t35;
	t1 = [-t28 * pkin(1) + t13 * r_i_i_C(1) - t39 * t9 - t31 * t10 + (r_i_i_C(2) * t20 + t40) * t36, -t31 * t11 + t39 * t12, -t12 * t15 + t16 * t37 + t42, t11, t42, 0; t30 * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t39 * t11 + t12 * t14 + t40 * t37, t39 * t10 - t31 * t9, -t10 * t15 - t16 * t36 + t43, t9, t43, 0; 0, (t39 * t27 + t31 * t29) * t25, -t15 * t38 + t26 * t16 + t41, -t25 * t29, t41, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (455->60), mult. (529->90), div. (0->0), fcn. (662->14), ass. (0->44)
	t62 = pkin(11) + r_i_i_C(3);
	t38 = sin(qJ(6));
	t41 = cos(qJ(6));
	t66 = t41 * r_i_i_C(1) - t38 * r_i_i_C(2) + pkin(5);
	t39 = sin(qJ(2));
	t40 = sin(qJ(1));
	t42 = cos(qJ(2));
	t43 = cos(qJ(1));
	t53 = cos(pkin(6));
	t50 = t43 * t53;
	t21 = t39 * t50 + t40 * t42;
	t36 = qJ(3) + pkin(12);
	t33 = qJ(5) + t36;
	t31 = sin(t33);
	t32 = cos(t33);
	t37 = sin(pkin(6));
	t54 = t37 * t43;
	t10 = t21 * t32 - t31 * t54;
	t20 = t40 * t39 - t42 * t50;
	t65 = t10 * t38 - t20 * t41;
	t64 = -t10 * t41 - t20 * t38;
	t28 = pkin(4) * cos(t36) + cos(qJ(3)) * pkin(3);
	t26 = pkin(2) + t28;
	t63 = t62 * t31 + t66 * t32 + t26;
	t57 = t37 * t39;
	t56 = t37 * t40;
	t55 = t37 * t42;
	t27 = pkin(4) * sin(t36) + sin(qJ(3)) * pkin(3);
	t52 = t37 * (pkin(8) + t27);
	t51 = t40 * t53;
	t35 = -pkin(10) - qJ(4) - pkin(9);
	t48 = t38 * r_i_i_C(1) + t41 * r_i_i_C(2) - t35;
	t9 = -t21 * t31 - t32 * t54;
	t47 = t62 * t10 + t66 * t9;
	t23 = -t39 * t51 + t43 * t42;
	t13 = t23 * t31 - t32 * t56;
	t14 = t23 * t32 + t31 * t56;
	t46 = -t66 * t13 + t62 * t14;
	t19 = t53 * t31 + t32 * t57;
	t45 = t62 * t19 + t66 * (-t31 * t57 + t53 * t32);
	t22 = t43 * t39 + t42 * t51;
	t2 = t14 * t41 + t22 * t38;
	t1 = -t14 * t38 + t22 * t41;
	t3 = [-t40 * pkin(1) - t10 * pkin(5) + t64 * r_i_i_C(1) + t65 * r_i_i_C(2) + t20 * t35 - t21 * t26 + t43 * t52 + t62 * t9, -t22 * t63 + t48 * t23, -t23 * t27 + t28 * t56 + t46, t22, t46, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t43 * pkin(1) + t14 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t62 * t13 - t22 * t35 + t23 * t26 + t40 * t52, -t20 * t63 + t48 * t21, -t21 * t27 - t28 * t54 + t47, t20, t47, -t65 * r_i_i_C(1) + t64 * r_i_i_C(2); 0, (t48 * t39 + t63 * t42) * t37, -t27 * t57 + t53 * t28 + t45, -t55, t45, (-t19 * t38 - t41 * t55) * r_i_i_C(1) + (-t19 * t41 + t38 * t55) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end