% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:44
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR9_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR9_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR9_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:44:09
	% EndTime: 2019-10-10 12:44:09
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:44:09
	% EndTime: 2019-10-10 12:44:09
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
	% StartTime: 2019-10-10 12:44:09
	% EndTime: 2019-10-10 12:44:09
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
	% StartTime: 2019-10-10 12:44:09
	% EndTime: 2019-10-10 12:44:10
	% DurationCPUTime: 0.15s
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
	% StartTime: 2019-10-10 12:44:09
	% EndTime: 2019-10-10 12:44:10
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (132->38), mult. (218->61), div. (0->0), fcn. (264->10), ass. (0->33)
	t20 = cos(pkin(6));
	t22 = sin(qJ(2));
	t26 = cos(qJ(1));
	t31 = t26 * t22;
	t23 = sin(qJ(1));
	t25 = cos(qJ(2));
	t32 = t23 * t25;
	t10 = t20 * t31 + t32;
	t18 = qJ(3) + qJ(4);
	t16 = sin(t18);
	t19 = sin(pkin(6));
	t34 = t19 * t26;
	t13 = t16 * t34;
	t17 = cos(t18);
	t40 = (-t10 * t16 - t17 * t34) * r_i_i_C(1) + (-t10 * t17 + t13) * r_i_i_C(2);
	t30 = t26 * t25;
	t33 = t23 * t22;
	t12 = -t20 * t33 + t30;
	t35 = t19 * t23;
	t5 = -t12 * t16 + t17 * t35;
	t6 = t12 * t17 + t16 * t35;
	t39 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t36 = t19 * t22;
	t38 = (-t16 * t36 + t20 * t17) * r_i_i_C(1) + (-t20 * t16 - t17 * t36) * r_i_i_C(2);
	t37 = r_i_i_C(3) + pkin(10) + pkin(9);
	t21 = sin(qJ(3));
	t29 = pkin(3) * t21 + pkin(8);
	t24 = cos(qJ(3));
	t15 = t24 * pkin(3) + pkin(2);
	t28 = t17 * r_i_i_C(1) - t16 * r_i_i_C(2) + t15;
	t11 = t20 * t32 + t31;
	t9 = -t20 * t30 + t33;
	t1 = [-t23 * pkin(1) + t13 * r_i_i_C(1) - t37 * t9 - t28 * t10 + (r_i_i_C(2) * t17 + t29) * t34, -t11 * t28 + t12 * t37, (-t12 * t21 + t24 * t35) * pkin(3) + t39, t39, 0, 0; t26 * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t11 * t37 + t12 * t15 + t29 * t35, t10 * t37 - t28 * t9, (-t10 * t21 - t24 * t34) * pkin(3) + t40, t40, 0, 0; 0, (t22 * t37 + t25 * t28) * t19, (t20 * t24 - t21 * t36) * pkin(3) + t38, t38, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:44:09
	% EndTime: 2019-10-10 12:44:10
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (290->51), mult. (462->80), div. (0->0), fcn. (578->12), ass. (0->38)
	t50 = r_i_i_C(3) + qJ(5);
	t30 = sin(pkin(12));
	t32 = cos(pkin(12));
	t57 = r_i_i_C(1) * t32 - r_i_i_C(2) * t30 + pkin(4);
	t36 = cos(qJ(3));
	t26 = t36 * pkin(3) + pkin(2);
	t29 = qJ(3) + qJ(4);
	t27 = sin(t29);
	t28 = cos(t29);
	t56 = t27 * t50 + t57 * t28 + t26;
	t31 = sin(pkin(6));
	t34 = sin(qJ(2));
	t53 = t31 * t34;
	t35 = sin(qJ(1));
	t52 = t31 * t35;
	t38 = cos(qJ(1));
	t51 = t31 * t38;
	t49 = cos(pkin(6));
	t37 = cos(qJ(2));
	t47 = t38 * t49;
	t19 = t34 * t47 + t35 * t37;
	t8 = t19 * t28 - t27 * t51;
	t48 = t35 * t49;
	t33 = sin(qJ(3));
	t46 = t31 * (pkin(3) * t33 + pkin(8));
	t39 = -pkin(10) - pkin(9);
	t44 = t30 * r_i_i_C(1) + t32 * r_i_i_C(2) - t39;
	t7 = t19 * t27 + t28 * t51;
	t43 = t50 * t8 - t57 * t7;
	t21 = -t34 * t48 + t38 * t37;
	t11 = t21 * t27 - t28 * t52;
	t12 = t21 * t28 + t27 * t52;
	t42 = -t57 * t11 + t50 * t12;
	t16 = t27 * t53 - t28 * t49;
	t41 = t50 * (t27 * t49 + t28 * t53) - t57 * t16;
	t20 = t38 * t34 + t37 * t48;
	t18 = t35 * t34 - t37 * t47;
	t1 = [(-t18 * t30 - t32 * t8) * r_i_i_C(1) + (-t18 * t32 + t30 * t8) * r_i_i_C(2) - t8 * pkin(4) - t19 * t26 + t18 * t39 - t35 * pkin(1) - t50 * t7 + t38 * t46, -t20 * t56 + t21 * t44, (-t21 * t33 + t36 * t52) * pkin(3) + t42, t42, t11, 0; (t12 * t32 + t20 * t30) * r_i_i_C(1) + (-t12 * t30 + t20 * t32) * r_i_i_C(2) + t12 * pkin(4) + t21 * t26 - t20 * t39 + t38 * pkin(1) + t35 * t46 + t50 * t11, -t18 * t56 + t19 * t44, (-t19 * t33 - t36 * t51) * pkin(3) + t43, t43, t7, 0; 0, (t44 * t34 + t56 * t37) * t31, (-t33 * t53 + t49 * t36) * pkin(3) + t41, t41, t16, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:44:09
	% EndTime: 2019-10-10 12:44:10
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (400->57), mult. (558->91), div. (0->0), fcn. (698->14), ass. (0->44)
	t64 = r_i_i_C(3) + pkin(11) + qJ(5);
	t31 = cos(pkin(12)) * pkin(5) + pkin(4);
	t37 = pkin(12) + qJ(6);
	t33 = sin(t37);
	t34 = cos(t37);
	t54 = t34 * r_i_i_C(1) - t33 * r_i_i_C(2) + t31;
	t45 = cos(qJ(3));
	t32 = t45 * pkin(3) + pkin(2);
	t38 = qJ(3) + qJ(4);
	t35 = sin(t38);
	t36 = cos(t38);
	t67 = t64 * t35 + t54 * t36 + t32;
	t40 = sin(pkin(6));
	t43 = sin(qJ(2));
	t63 = t40 * t43;
	t44 = sin(qJ(1));
	t62 = t40 * t44;
	t46 = cos(qJ(2));
	t61 = t40 * t46;
	t47 = cos(qJ(1));
	t60 = t40 * t47;
	t59 = cos(pkin(6));
	t58 = sin(pkin(12)) * pkin(5) + pkin(10) + pkin(9);
	t56 = t47 * t59;
	t24 = t43 * t56 + t44 * t46;
	t12 = t24 * t36 - t35 * t60;
	t57 = t44 * t59;
	t42 = sin(qJ(3));
	t55 = t40 * (pkin(3) * t42 + pkin(8));
	t11 = t24 * t35 + t36 * t60;
	t53 = -t54 * t11 + t64 * t12;
	t26 = -t43 * t57 + t47 * t46;
	t15 = t26 * t35 - t36 * t62;
	t16 = t26 * t36 + t35 * t62;
	t52 = -t54 * t15 + t64 * t16;
	t21 = t35 * t63 - t59 * t36;
	t22 = t59 * t35 + t36 * t63;
	t51 = -t54 * t21 + t64 * t22;
	t50 = t33 * r_i_i_C(1) + t34 * r_i_i_C(2) + t58;
	t25 = t47 * t43 + t46 * t57;
	t23 = t44 * t43 - t46 * t56;
	t2 = t16 * t34 + t25 * t33;
	t1 = -t16 * t33 + t25 * t34;
	t3 = [-t44 * pkin(1) - t64 * t11 - t54 * t12 - t50 * t23 - t24 * t32 + t47 * t55, -t25 * t67 + t50 * t26, (-t26 * t42 + t45 * t62) * pkin(3) + t52, t52, t15, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t47 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t64 * t15 + t16 * t31 + t58 * t25 + t26 * t32 + t44 * t55, -t23 * t67 + t50 * t24, (-t24 * t42 - t45 * t60) * pkin(3) + t53, t53, t11, (-t12 * t33 + t23 * t34) * r_i_i_C(1) + (-t12 * t34 - t23 * t33) * r_i_i_C(2); 0, (t50 * t43 + t67 * t46) * t40, (-t42 * t63 + t59 * t45) * pkin(3) + t51, t51, t21, (-t22 * t33 - t34 * t61) * r_i_i_C(1) + (-t22 * t34 + t33 * t61) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end