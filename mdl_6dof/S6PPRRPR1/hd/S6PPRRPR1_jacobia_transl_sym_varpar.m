% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRPR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRPR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobia_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t1 = sin(pkin(6));
	t2 = [0, sin(pkin(11)) * t1, 0, 0, 0, 0; 0, -cos(pkin(11)) * t1, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (24->18), mult. (70->40), div. (0->0), fcn. (93->10), ass. (0->20)
	t6 = sin(pkin(11));
	t8 = sin(pkin(6));
	t21 = t6 * t8;
	t10 = cos(pkin(11));
	t20 = t10 * t8;
	t11 = cos(pkin(7));
	t9 = cos(pkin(12));
	t19 = t11 * t9;
	t12 = cos(pkin(6));
	t18 = t12 * t6;
	t17 = t10 * t12;
	t5 = sin(pkin(12));
	t7 = sin(pkin(7));
	t16 = t11 * (-t10 * t5 - t9 * t18) + t7 * t21;
	t15 = -(t9 * t17 - t5 * t6) * t11 + t7 * t20;
	t14 = cos(qJ(3));
	t13 = sin(qJ(3));
	t4 = t10 * t9 - t5 * t18;
	t2 = t5 * t17 + t6 * t9;
	t1 = [0, t21, (-t4 * t13 + t16 * t14) * r_i_i_C(1) + (-t16 * t13 - t14 * t4) * r_i_i_C(2), 0, 0, 0; 0, -t20, (-t2 * t13 - t15 * t14) * r_i_i_C(1) + (t15 * t13 - t14 * t2) * r_i_i_C(2), 0, 0, 0; 1, t12, (r_i_i_C(1) * t14 - r_i_i_C(2) * t13) * t7 * t12 + ((-t13 * t5 + t14 * t19) * r_i_i_C(1) + (-t13 * t19 - t14 * t5) * r_i_i_C(2)) * t8, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (102->33), mult. (288->67), div. (0->0), fcn. (377->12), ass. (0->34)
	t36 = pkin(9) + r_i_i_C(3);
	t15 = sin(pkin(11));
	t17 = sin(pkin(6));
	t35 = t15 * t17;
	t21 = cos(pkin(6));
	t34 = t15 * t21;
	t16 = sin(pkin(7));
	t33 = t16 * t17;
	t32 = t16 * t21;
	t18 = cos(pkin(12));
	t20 = cos(pkin(7));
	t31 = t18 * t20;
	t19 = cos(pkin(11));
	t30 = t19 * t17;
	t29 = t19 * t21;
	t22 = sin(qJ(4));
	t24 = cos(qJ(4));
	t28 = t24 * r_i_i_C(1) - t22 * r_i_i_C(2) + pkin(3);
	t14 = sin(pkin(12));
	t10 = -t15 * t14 + t18 * t29;
	t27 = t10 * t20 - t16 * t30;
	t12 = -t19 * t14 - t18 * t34;
	t26 = t12 * t20 + t15 * t33;
	t25 = cos(qJ(3));
	t23 = sin(qJ(3));
	t13 = -t14 * t34 + t19 * t18;
	t11 = t14 * t29 + t15 * t18;
	t9 = -t18 * t33 + t21 * t20;
	t8 = -t12 * t16 + t20 * t35;
	t7 = -t10 * t16 - t20 * t30;
	t6 = t23 * t32 + (t14 * t25 + t23 * t31) * t17;
	t4 = t13 * t25 + t26 * t23;
	t2 = t11 * t25 + t27 * t23;
	t1 = [0, t35, t36 * t4 + t28 * (-t13 * t23 + t26 * t25), (-t4 * t22 + t8 * t24) * r_i_i_C(1) + (-t8 * t22 - t4 * t24) * r_i_i_C(2), 0, 0; 0, -t30, t36 * t2 + t28 * (-t11 * t23 + t27 * t25), (-t2 * t22 + t7 * t24) * r_i_i_C(1) + (-t2 * t24 - t7 * t22) * r_i_i_C(2), 0, 0; 1, t21, t36 * t6 + t28 * (t25 * t32 + (-t14 * t23 + t25 * t31) * t17), (-t6 * t22 + t9 * t24) * r_i_i_C(1) + (-t9 * t22 - t6 * t24) * r_i_i_C(2), 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (235->37), mult. (662->71), div. (0->0), fcn. (872->14), ass. (0->41)
	t22 = sin(pkin(11));
	t24 = sin(pkin(6));
	t49 = t22 * t24;
	t29 = cos(pkin(6));
	t48 = t22 * t29;
	t23 = sin(pkin(7));
	t47 = t23 * t24;
	t46 = t23 * t29;
	t26 = cos(pkin(12));
	t28 = cos(pkin(7));
	t45 = t26 * t28;
	t27 = cos(pkin(11));
	t44 = t27 * t24;
	t43 = t27 * t29;
	t42 = r_i_i_C(3) + qJ(5);
	t20 = sin(pkin(13));
	t25 = cos(pkin(13));
	t41 = r_i_i_C(1) * t25 - r_i_i_C(2) * t20 + pkin(4);
	t40 = t20 * r_i_i_C(1) + t25 * r_i_i_C(2) + pkin(9);
	t21 = sin(pkin(12));
	t16 = -t22 * t21 + t26 * t43;
	t39 = -t16 * t23 - t28 * t44;
	t38 = t16 * t28 - t23 * t44;
	t18 = -t27 * t21 - t26 * t48;
	t37 = -t18 * t23 + t28 * t49;
	t36 = t18 * t28 + t22 * t47;
	t35 = -t26 * t47 + t29 * t28;
	t30 = sin(qJ(4));
	t32 = cos(qJ(4));
	t34 = t42 * t30 + t41 * t32 + pkin(3);
	t33 = cos(qJ(3));
	t31 = sin(qJ(3));
	t19 = -t21 * t48 + t27 * t26;
	t17 = t21 * t43 + t22 * t26;
	t14 = t31 * t46 + (t21 * t33 + t31 * t45) * t24;
	t9 = t14 * t30 - t35 * t32;
	t8 = t19 * t33 + t36 * t31;
	t6 = t17 * t33 + t38 * t31;
	t3 = t8 * t30 - t37 * t32;
	t1 = t6 * t30 - t39 * t32;
	t2 = [0, t49, t40 * t8 + t34 * (-t19 * t31 + t36 * t33), t42 * (t37 * t30 + t8 * t32) - t41 * t3, t3, 0; 0, -t44, t40 * t6 + t34 * (-t17 * t31 + t38 * t33), t42 * (t39 * t30 + t6 * t32) - t41 * t1, t1, 0; 1, t29, t40 * t14 + t34 * (t33 * t46 + (-t21 * t31 + t33 * t45) * t24), t42 * (t14 * t32 + t35 * t30) - t41 * t9, t9, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (355->48), mult. (889->87), div. (0->0), fcn. (1171->16), ass. (0->50)
	t54 = sin(pkin(12));
	t55 = sin(pkin(11));
	t45 = t55 * t54;
	t58 = cos(pkin(12));
	t59 = cos(pkin(11));
	t52 = t59 * t58;
	t61 = cos(pkin(6));
	t37 = -t61 * t52 + t45;
	t56 = sin(pkin(7));
	t57 = sin(pkin(6));
	t49 = t57 * t56;
	t60 = cos(pkin(7));
	t66 = t37 * t60 + t59 * t49;
	t47 = t55 * t58;
	t50 = t59 * t54;
	t38 = t61 * t47 + t50;
	t46 = t55 * t57;
	t65 = t38 * t60 - t56 * t46;
	t64 = t58 * t60 * t57 + t61 * t56;
	t63 = r_i_i_C(3) + pkin(10) + qJ(5);
	t62 = cos(qJ(3));
	t51 = t59 * t57;
	t48 = t57 * t54;
	t26 = pkin(13) + qJ(6);
	t24 = sin(t26);
	t25 = cos(t26);
	t44 = -t25 * r_i_i_C(1) + t24 * r_i_i_C(2) - cos(pkin(13)) * pkin(5) - pkin(4);
	t43 = sin(pkin(13)) * pkin(5) + t24 * r_i_i_C(1) + t25 * r_i_i_C(2) + pkin(9);
	t29 = sin(qJ(4));
	t31 = cos(qJ(4));
	t39 = -t63 * t29 + t44 * t31 - pkin(3);
	t36 = -t58 * t49 + t61 * t60;
	t33 = t38 * t56 + t60 * t46;
	t32 = t37 * t56 - t60 * t51;
	t30 = sin(qJ(3));
	t19 = -t61 * t45 + t52;
	t18 = t61 * t50 + t47;
	t14 = t64 * t30 + t62 * t48;
	t13 = t30 * t48 - t64 * t62;
	t10 = t14 * t31 + t36 * t29;
	t9 = t14 * t29 - t36 * t31;
	t8 = t19 * t62 - t65 * t30;
	t7 = t19 * t30 + t65 * t62;
	t6 = t18 * t62 - t66 * t30;
	t5 = t18 * t30 + t66 * t62;
	t4 = t33 * t29 + t8 * t31;
	t3 = t8 * t29 - t33 * t31;
	t2 = t32 * t29 + t6 * t31;
	t1 = t6 * t29 - t32 * t31;
	t11 = [0, t46, t39 * t7 + t43 * t8, t44 * t3 + t63 * t4, t3, (-t4 * t24 + t7 * t25) * r_i_i_C(1) + (-t7 * t24 - t4 * t25) * r_i_i_C(2); 0, -t51, t39 * t5 + t43 * t6, t44 * t1 + t63 * t2, t1, (-t2 * t24 + t5 * t25) * r_i_i_C(1) + (-t2 * t25 - t5 * t24) * r_i_i_C(2); 1, t61, t39 * t13 + t43 * t14, t63 * t10 + t44 * t9, t9, (-t10 * t24 + t13 * t25) * r_i_i_C(1) + (-t10 * t25 - t13 * t24) * r_i_i_C(2);];
	Ja_transl = t11;
else
	Ja_transl=NaN(3,6);
end