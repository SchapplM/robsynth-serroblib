% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:12
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRPR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRPR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t1 = sin(pkin(6));
	t2 = [0, sin(pkin(11)) * t1, 0, 0, 0, 0; 0, -cos(pkin(11)) * t1, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.13s
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
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.18s
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
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (188->35), mult. (524->67), div. (0->0), fcn. (690->12), ass. (0->39)
	t44 = pkin(4) - r_i_i_C(2);
	t43 = pkin(9) + r_i_i_C(1);
	t21 = sin(pkin(11));
	t23 = sin(pkin(6));
	t42 = t21 * t23;
	t27 = cos(pkin(6));
	t41 = t21 * t27;
	t22 = sin(pkin(7));
	t40 = t22 * t23;
	t39 = t22 * t27;
	t24 = cos(pkin(12));
	t26 = cos(pkin(7));
	t38 = t24 * t26;
	t25 = cos(pkin(11));
	t37 = t25 * t23;
	t36 = t25 * t27;
	t35 = r_i_i_C(3) + qJ(5);
	t20 = sin(pkin(12));
	t16 = -t21 * t20 + t24 * t36;
	t34 = t16 * t26 - t22 * t37;
	t18 = -t25 * t20 - t24 * t41;
	t33 = t18 * t26 + t21 * t40;
	t28 = sin(qJ(4));
	t30 = cos(qJ(4));
	t32 = t35 * t28 + t44 * t30 + pkin(3);
	t31 = cos(qJ(3));
	t29 = sin(qJ(3));
	t19 = -t20 * t41 + t25 * t24;
	t17 = t20 * t36 + t21 * t24;
	t15 = -t24 * t40 + t27 * t26;
	t14 = -t18 * t22 + t26 * t42;
	t13 = -t16 * t22 - t26 * t37;
	t12 = t29 * t39 + (t20 * t31 + t29 * t38) * t23;
	t9 = t12 * t28 - t15 * t30;
	t8 = t19 * t31 + t33 * t29;
	t6 = t17 * t31 + t34 * t29;
	t3 = -t14 * t30 + t8 * t28;
	t1 = -t13 * t30 + t6 * t28;
	t2 = [0, t42, t43 * t8 + t32 * (-t19 * t29 + t33 * t31), t35 * (t14 * t28 + t8 * t30) - t44 * t3, t3, 0; 0, -t37, t43 * t6 + t32 * (-t17 * t29 + t34 * t31), t35 * (t13 * t28 + t6 * t30) - t44 * t1, t1, 0; 1, t27, t43 * t12 + t32 * (t31 * t39 + (-t20 * t29 + t31 * t38) * t23), -t44 * t9 + t35 * (t12 * t30 + t15 * t28), t9, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (347->46), mult. (972->85), div. (0->0), fcn. (1282->14), ass. (0->46)
	t51 = sin(pkin(12));
	t52 = sin(pkin(11));
	t41 = t52 * t51;
	t55 = cos(pkin(12));
	t56 = cos(pkin(11));
	t49 = t56 * t55;
	t58 = cos(pkin(6));
	t33 = -t58 * t49 + t41;
	t53 = sin(pkin(7));
	t54 = sin(pkin(6));
	t46 = t54 * t53;
	t57 = cos(pkin(7));
	t62 = t33 * t57 + t56 * t46;
	t43 = t52 * t55;
	t47 = t56 * t51;
	t34 = t58 * t43 + t47;
	t42 = t52 * t54;
	t61 = t34 * t57 - t53 * t42;
	t60 = t55 * t57 * t54 + t53 * t58;
	t59 = cos(qJ(3));
	t50 = -pkin(4) - pkin(10) - r_i_i_C(3);
	t48 = t56 * t54;
	t45 = t54 * t51;
	t23 = sin(qJ(6));
	t26 = cos(qJ(6));
	t40 = t23 * r_i_i_C(1) + t26 * r_i_i_C(2) + qJ(5);
	t39 = t26 * r_i_i_C(1) - t23 * r_i_i_C(2) + pkin(5) + pkin(9);
	t24 = sin(qJ(4));
	t27 = cos(qJ(4));
	t35 = -t40 * t24 + t50 * t27 - pkin(3);
	t32 = -t55 * t46 + t58 * t57;
	t29 = t34 * t53 + t57 * t42;
	t28 = t33 * t53 - t57 * t48;
	t25 = sin(qJ(3));
	t19 = -t58 * t41 + t49;
	t18 = t58 * t47 + t43;
	t14 = t60 * t25 + t59 * t45;
	t13 = t25 * t45 - t60 * t59;
	t9 = t14 * t24 - t32 * t27;
	t8 = t19 * t59 - t61 * t25;
	t7 = t19 * t25 + t61 * t59;
	t6 = t18 * t59 - t62 * t25;
	t5 = t18 * t25 + t62 * t59;
	t3 = t8 * t24 - t29 * t27;
	t1 = t6 * t24 - t28 * t27;
	t2 = [0, t42, t35 * t7 + t39 * t8, t40 * (t29 * t24 + t8 * t27) + t50 * t3, t3, (-t7 * t23 + t3 * t26) * r_i_i_C(1) + (-t3 * t23 - t7 * t26) * r_i_i_C(2); 0, -t48, t35 * t5 + t39 * t6, t40 * (t28 * t24 + t6 * t27) + t50 * t1, t1, (t1 * t26 - t5 * t23) * r_i_i_C(1) + (-t1 * t23 - t5 * t26) * r_i_i_C(2); 1, t58, t35 * t13 + t39 * t14, t50 * t9 + t40 * (t14 * t27 + t32 * t24), t9, (-t13 * t23 + t9 * t26) * r_i_i_C(1) + (-t13 * t26 - t9 * t23) * r_i_i_C(2);];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,6);
end