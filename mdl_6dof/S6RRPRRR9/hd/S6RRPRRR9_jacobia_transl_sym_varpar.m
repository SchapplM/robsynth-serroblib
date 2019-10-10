% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR9
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
% Datum: 2019-10-10 11:03
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR9_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR9_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR9_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
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
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
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
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
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
	t6 = sin(pkin(12));
	t8 = cos(pkin(12));
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
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
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
	t22 = sin(pkin(12)) * pkin(3) + pkin(8);
	t12 = pkin(12) + qJ(4);
	t10 = sin(t12);
	t11 = cos(t12);
	t9 = cos(pkin(12)) * pkin(3) + pkin(2);
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
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (198->42), mult. (232->64), div. (0->0), fcn. (283->12), ass. (0->34)
	t25 = cos(pkin(6));
	t26 = sin(qJ(2));
	t29 = cos(qJ(1));
	t32 = t29 * t26;
	t27 = sin(qJ(1));
	t28 = cos(qJ(2));
	t33 = t27 * t28;
	t10 = t25 * t32 + t33;
	t23 = pkin(12) + qJ(4);
	t21 = qJ(5) + t23;
	t17 = sin(t21);
	t24 = sin(pkin(6));
	t35 = t24 * t29;
	t14 = t17 * t35;
	t18 = cos(t21);
	t42 = (-t10 * t17 - t18 * t35) * r_i_i_C(1) + (-t10 * t18 + t14) * r_i_i_C(2);
	t31 = t29 * t28;
	t34 = t27 * t26;
	t12 = -t25 * t34 + t31;
	t36 = t24 * t27;
	t5 = -t12 * t17 + t18 * t36;
	t6 = t12 * t18 + t17 * t36;
	t41 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t37 = t24 * t26;
	t40 = (-t17 * t37 + t25 * t18) * r_i_i_C(1) + (-t25 * t17 - t18 * t37) * r_i_i_C(2);
	t19 = sin(t23);
	t39 = pkin(8) + pkin(4) * t19 + sin(pkin(12)) * pkin(3);
	t38 = r_i_i_C(3) + pkin(10) + pkin(9) + qJ(3);
	t20 = cos(t23);
	t13 = pkin(4) * t20 + cos(pkin(12)) * pkin(3) + pkin(2);
	t30 = t18 * r_i_i_C(1) - t17 * r_i_i_C(2) + t13;
	t11 = t25 * t33 + t32;
	t9 = -t25 * t31 + t34;
	t1 = [-t27 * pkin(1) + t14 * r_i_i_C(1) - t38 * t9 - t30 * t10 + (r_i_i_C(2) * t18 + t39) * t35, -t30 * t11 + t38 * t12, t11, (-t12 * t19 + t20 * t36) * pkin(4) + t41, t41, 0; t29 * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t38 * t11 + t12 * t13 + t39 * t36, t38 * t10 - t30 * t9, t9, (-t10 * t19 - t20 * t35) * pkin(4) + t42, t42, 0; 0, (t38 * t26 + t30 * t28) * t24, -t24 * t28, (-t19 * t37 + t20 * t25) * pkin(4) + t40, t40, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (449->60), mult. (523->93), div. (0->0), fcn. (656->14), ass. (0->44)
	t61 = pkin(11) + r_i_i_C(3);
	t37 = sin(qJ(6));
	t40 = cos(qJ(6));
	t65 = t40 * r_i_i_C(1) - t37 * r_i_i_C(2) + pkin(5);
	t38 = sin(qJ(2));
	t39 = sin(qJ(1));
	t41 = cos(qJ(2));
	t42 = cos(qJ(1));
	t52 = cos(pkin(6));
	t49 = t42 * t52;
	t21 = t38 * t49 + t39 * t41;
	t35 = pkin(12) + qJ(4);
	t33 = qJ(5) + t35;
	t29 = sin(t33);
	t30 = cos(t33);
	t36 = sin(pkin(6));
	t53 = t36 * t42;
	t10 = t21 * t30 - t29 * t53;
	t20 = t39 * t38 - t41 * t49;
	t64 = t10 * t37 - t20 * t40;
	t63 = -t10 * t40 - t20 * t37;
	t32 = cos(t35);
	t24 = pkin(4) * t32 + cos(pkin(12)) * pkin(3) + pkin(2);
	t62 = t61 * t29 + t65 * t30 + t24;
	t56 = t36 * t38;
	t55 = t36 * t39;
	t54 = t36 * t41;
	t31 = sin(t35);
	t51 = t36 * (pkin(8) + pkin(4) * t31 + sin(pkin(12)) * pkin(3));
	t50 = t39 * t52;
	t34 = -pkin(10) - pkin(9) - qJ(3);
	t47 = t37 * r_i_i_C(1) + t40 * r_i_i_C(2) - t34;
	t9 = -t21 * t29 - t30 * t53;
	t46 = t61 * t10 + t65 * t9;
	t23 = -t38 * t50 + t42 * t41;
	t13 = t23 * t29 - t30 * t55;
	t14 = t23 * t30 + t29 * t55;
	t45 = -t65 * t13 + t61 * t14;
	t19 = t52 * t29 + t30 * t56;
	t44 = t61 * t19 + t65 * (-t29 * t56 + t52 * t30);
	t22 = t42 * t38 + t41 * t50;
	t2 = t14 * t40 + t22 * t37;
	t1 = -t14 * t37 + t22 * t40;
	t3 = [-t39 * pkin(1) - t10 * pkin(5) + t63 * r_i_i_C(1) + t64 * r_i_i_C(2) + t20 * t34 - t21 * t24 + t42 * t51 + t61 * t9, -t22 * t62 + t47 * t23, t22, (-t23 * t31 + t32 * t55) * pkin(4) + t45, t45, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t42 * pkin(1) + t14 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t61 * t13 - t22 * t34 + t23 * t24 + t39 * t51, -t20 * t62 + t47 * t21, t20, (-t21 * t31 - t32 * t53) * pkin(4) + t46, t46, -t64 * r_i_i_C(1) + t63 * r_i_i_C(2); 0, (t47 * t38 + t62 * t41) * t36, -t54, (-t31 * t56 + t52 * t32) * pkin(4) + t44, t44, (-t19 * t37 - t40 * t54) * r_i_i_C(1) + (-t19 * t40 + t37 * t54) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end