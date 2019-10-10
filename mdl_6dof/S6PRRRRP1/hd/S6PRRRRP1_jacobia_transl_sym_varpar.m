% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:02
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRP1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:02:05
	% EndTime: 2019-10-09 23:02:05
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:02:05
	% EndTime: 2019-10-09 23:02:05
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:02:05
	% EndTime: 2019-10-09 23:02:05
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(6));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(11));
	t1 = sin(pkin(11));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:02:05
	% EndTime: 2019-10-09 23:02:05
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (35->19), mult. (92->37), div. (0->0), fcn. (112->8), ass. (0->19)
	t6 = sin(pkin(6));
	t9 = sin(qJ(3));
	t20 = t6 * t9;
	t19 = pkin(8) + r_i_i_C(3);
	t11 = cos(qJ(3));
	t18 = t11 * t6;
	t10 = sin(qJ(2));
	t5 = sin(pkin(11));
	t17 = t5 * t10;
	t12 = cos(qJ(2));
	t16 = t5 * t12;
	t7 = cos(pkin(11));
	t15 = t7 * t10;
	t14 = t7 * t12;
	t13 = t11 * r_i_i_C(1) - t9 * r_i_i_C(2) + pkin(2);
	t8 = cos(pkin(6));
	t4 = -t8 * t17 + t14;
	t2 = t8 * t15 + t16;
	t1 = [0, t19 * t4 + t13 * (-t8 * t16 - t15), (t5 * t18 - t4 * t9) * r_i_i_C(1) + (-t11 * t4 - t5 * t20) * r_i_i_C(2), 0, 0, 0; 0, t19 * t2 + t13 * (t8 * t14 - t17), (-t7 * t18 - t2 * t9) * r_i_i_C(1) + (-t11 * t2 + t7 * t20) * r_i_i_C(2), 0, 0, 0; 1, (t19 * t10 + t13 * t12) * t6, (-t10 * t20 + t11 * t8) * r_i_i_C(1) + (-t10 * t18 - t8 * t9) * r_i_i_C(2), 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:02:05
	% EndTime: 2019-10-09 23:02:05
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (92->27), mult. (152->51), div. (0->0), fcn. (184->10), ass. (0->25)
	t14 = qJ(3) + qJ(4);
	t12 = sin(t14);
	t13 = cos(t14);
	t16 = sin(pkin(6));
	t17 = cos(pkin(11));
	t29 = t16 * t17;
	t15 = sin(pkin(11));
	t22 = cos(qJ(2));
	t18 = cos(pkin(6));
	t20 = sin(qJ(2));
	t26 = t18 * t20;
	t8 = t15 * t22 + t17 * t26;
	t34 = (-t12 * t8 - t13 * t29) * r_i_i_C(1) + (t12 * t29 - t13 * t8) * r_i_i_C(2);
	t10 = -t15 * t26 + t17 * t22;
	t30 = t15 * t16;
	t33 = (-t10 * t12 + t13 * t30) * r_i_i_C(1) + (-t10 * t13 - t12 * t30) * r_i_i_C(2);
	t28 = t16 * t20;
	t32 = (-t12 * t28 + t13 * t18) * r_i_i_C(1) + (-t12 * t18 - t13 * t28) * r_i_i_C(2);
	t31 = r_i_i_C(3) + pkin(9) + pkin(8);
	t21 = cos(qJ(3));
	t27 = t16 * t21;
	t25 = t18 * t22;
	t24 = pkin(3) * t21 + r_i_i_C(1) * t13 - r_i_i_C(2) * t12 + pkin(2);
	t19 = sin(qJ(3));
	t1 = [0, t31 * t10 + t24 * (-t15 * t25 - t17 * t20), (-t10 * t19 + t15 * t27) * pkin(3) + t33, t33, 0, 0; 0, t31 * t8 + t24 * (-t15 * t20 + t17 * t25), (-t17 * t27 - t19 * t8) * pkin(3) + t34, t34, 0, 0; 1, (t31 * t20 + t24 * t22) * t16, (t18 * t21 - t19 * t28) * pkin(3) + t32, t32, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:02:05
	% EndTime: 2019-10-09 23:02:05
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (240->41), mult. (389->71), div. (0->0), fcn. (485->12), ass. (0->36)
	t56 = pkin(10) + r_i_i_C(3);
	t29 = sin(qJ(5));
	t32 = cos(qJ(5));
	t55 = t32 * r_i_i_C(1) - t29 * r_i_i_C(2) + pkin(4);
	t25 = qJ(3) + qJ(4);
	t23 = sin(t25);
	t24 = cos(t25);
	t33 = cos(qJ(3));
	t54 = t33 * pkin(3) + t56 * t23 + t55 * t24 + pkin(2);
	t26 = sin(pkin(11));
	t27 = sin(pkin(6));
	t50 = t26 * t27;
	t31 = sin(qJ(2));
	t49 = t26 * t31;
	t34 = cos(qJ(2));
	t48 = t26 * t34;
	t47 = t27 * t31;
	t46 = t27 * t34;
	t45 = cos(pkin(11));
	t44 = t27 * t45;
	t43 = t45 * t31;
	t42 = t45 * t34;
	t40 = t29 * r_i_i_C(1) + t32 * r_i_i_C(2) + pkin(8) + pkin(9);
	t28 = cos(pkin(6));
	t17 = t28 * t43 + t48;
	t8 = t17 * t24 - t23 * t44;
	t39 = t56 * t8 + t55 * (-t17 * t23 - t24 * t44);
	t19 = -t28 * t49 + t42;
	t10 = t19 * t24 + t23 * t50;
	t38 = t56 * t10 + t55 * (-t19 * t23 + t24 * t50);
	t15 = t28 * t23 + t24 * t47;
	t37 = t56 * t15 + t55 * (-t23 * t47 + t28 * t24);
	t30 = sin(qJ(3));
	t18 = t28 * t48 + t43;
	t16 = -t28 * t42 + t49;
	t1 = [0, -t18 * t54 + t40 * t19, (-t19 * t30 + t33 * t50) * pkin(3) + t38, t38, (-t10 * t29 + t18 * t32) * r_i_i_C(1) + (-t10 * t32 - t18 * t29) * r_i_i_C(2), 0; 0, -t16 * t54 + t40 * t17, (-t17 * t30 - t33 * t44) * pkin(3) + t39, t39, (t16 * t32 - t8 * t29) * r_i_i_C(1) + (-t16 * t29 - t8 * t32) * r_i_i_C(2), 0; 1, (t40 * t31 + t34 * t54) * t27, (t28 * t33 - t30 * t47) * pkin(3) + t37, t37, (-t15 * t29 - t32 * t46) * r_i_i_C(1) + (-t15 * t32 + t29 * t46) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:02:05
	% EndTime: 2019-10-09 23:02:05
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (293->42), mult. (460->74), div. (0->0), fcn. (573->12), ass. (0->39)
	t63 = pkin(5) + r_i_i_C(1);
	t65 = r_i_i_C(3) + qJ(6) + pkin(10);
	t35 = sin(qJ(5));
	t38 = cos(qJ(5));
	t64 = -t35 * r_i_i_C(2) + t63 * t38 + pkin(4);
	t30 = qJ(3) + qJ(4);
	t28 = sin(t30);
	t29 = cos(t30);
	t39 = cos(qJ(3));
	t62 = t39 * pkin(3) + t65 * t28 + t64 * t29 + pkin(2);
	t31 = sin(pkin(11));
	t32 = sin(pkin(6));
	t58 = t31 * t32;
	t33 = cos(pkin(11));
	t57 = t32 * t33;
	t37 = sin(qJ(2));
	t56 = t32 * t37;
	t55 = t32 * t39;
	t40 = cos(qJ(2));
	t54 = t32 * t40;
	t53 = cos(pkin(6));
	t52 = t37 * t53;
	t51 = t40 * t53;
	t20 = t31 * t40 + t33 * t52;
	t10 = t20 * t29 - t28 * t57;
	t9 = t20 * t28 + t29 * t57;
	t46 = t65 * t10 - t64 * t9;
	t22 = -t31 * t52 + t33 * t40;
	t11 = t22 * t28 - t29 * t58;
	t12 = t22 * t29 + t28 * t58;
	t45 = -t64 * t11 + t65 * t12;
	t44 = t38 * r_i_i_C(2) + t63 * t35 + pkin(8) + pkin(9);
	t17 = t28 * t56 - t53 * t29;
	t18 = t53 * t28 + t29 * t56;
	t43 = -t64 * t17 + t65 * t18;
	t36 = sin(qJ(3));
	t21 = t31 * t51 + t33 * t37;
	t19 = t31 * t37 - t33 * t51;
	t1 = [0, -t21 * t62 + t44 * t22, (-t22 * t36 + t31 * t55) * pkin(3) + t45, t45, (-t12 * t38 - t21 * t35) * r_i_i_C(2) + t63 * (-t12 * t35 + t21 * t38), t11; 0, -t19 * t62 + t44 * t20, (-t20 * t36 - t33 * t55) * pkin(3) + t46, t46, (-t10 * t38 - t19 * t35) * r_i_i_C(2) + t63 * (-t10 * t35 + t19 * t38), t9; 1, (t44 * t37 + t62 * t40) * t32, (-t36 * t56 + t53 * t39) * pkin(3) + t43, t43, (-t18 * t38 + t35 * t54) * r_i_i_C(2) + t63 * (-t18 * t35 - t38 * t54), t17;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end