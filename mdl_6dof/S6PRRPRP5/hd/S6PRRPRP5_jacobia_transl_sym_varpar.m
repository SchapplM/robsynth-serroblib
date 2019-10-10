% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRP5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:23
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRP5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRP5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:23:43
	% EndTime: 2019-10-09 22:23:43
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:23:43
	% EndTime: 2019-10-09 22:23:43
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:23:43
	% EndTime: 2019-10-09 22:23:43
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(6));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(10));
	t1 = sin(pkin(10));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:23:43
	% EndTime: 2019-10-09 22:23:43
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (35->19), mult. (92->37), div. (0->0), fcn. (112->8), ass. (0->19)
	t6 = sin(pkin(6));
	t9 = sin(qJ(3));
	t20 = t6 * t9;
	t19 = pkin(8) + r_i_i_C(3);
	t11 = cos(qJ(3));
	t18 = t11 * t6;
	t10 = sin(qJ(2));
	t5 = sin(pkin(10));
	t17 = t5 * t10;
	t12 = cos(qJ(2));
	t16 = t5 * t12;
	t7 = cos(pkin(10));
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
	% StartTime: 2019-10-09 22:23:43
	% EndTime: 2019-10-09 22:23:43
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (66->21), mult. (169->39), div. (0->0), fcn. (211->8), ass. (0->22)
	t26 = pkin(3) - r_i_i_C(2);
	t25 = pkin(8) + r_i_i_C(1);
	t12 = sin(pkin(6));
	t15 = sin(qJ(3));
	t24 = t12 * t15;
	t17 = cos(qJ(3));
	t23 = t12 * t17;
	t14 = cos(pkin(6));
	t16 = sin(qJ(2));
	t22 = t14 * t16;
	t18 = cos(qJ(2));
	t21 = t14 * t18;
	t20 = r_i_i_C(3) + qJ(4);
	t19 = t20 * t15 + t26 * t17 + pkin(2);
	t13 = cos(pkin(10));
	t11 = sin(pkin(10));
	t9 = -t14 * t17 + t16 * t24;
	t8 = -t11 * t22 + t13 * t18;
	t6 = t11 * t18 + t13 * t22;
	t3 = -t11 * t23 + t15 * t8;
	t1 = t13 * t23 + t15 * t6;
	t2 = [0, t25 * t8 + t19 * (-t11 * t21 - t13 * t16), t20 * (t11 * t24 + t17 * t8) - t26 * t3, t3, 0, 0; 0, t25 * t6 + t19 * (-t11 * t16 + t13 * t21), t20 * (-t13 * t24 + t17 * t6) - t26 * t1, t1, 0, 0; 1, (t16 * t25 + t18 * t19) * t12, -t26 * t9 + t20 * (t14 * t15 + t16 * t23), t9, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:23:43
	% EndTime: 2019-10-09 22:23:43
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (125->33), mult. (325->63), div. (0->0), fcn. (411->10), ass. (0->28)
	t17 = sin(qJ(3));
	t20 = cos(qJ(3));
	t16 = sin(qJ(5));
	t19 = cos(qJ(5));
	t24 = t16 * r_i_i_C(1) + t19 * r_i_i_C(2) + qJ(4);
	t28 = pkin(3) + pkin(9) + r_i_i_C(3);
	t34 = t24 * t17 + t28 * t20 + pkin(2);
	t15 = sin(pkin(6));
	t33 = t15 * t17;
	t32 = t15 * t20;
	t21 = cos(qJ(2));
	t31 = t15 * t21;
	t30 = cos(pkin(6));
	t29 = cos(pkin(10));
	t14 = sin(pkin(10));
	t27 = t14 * t30;
	t26 = t15 * t29;
	t25 = t30 * t29;
	t23 = t19 * r_i_i_C(1) - t16 * r_i_i_C(2) + pkin(4) + pkin(8);
	t18 = sin(qJ(2));
	t9 = t18 * t33 - t30 * t20;
	t8 = -t18 * t27 + t29 * t21;
	t7 = t29 * t18 + t21 * t27;
	t6 = t14 * t21 + t18 * t25;
	t5 = t14 * t18 - t21 * t25;
	t3 = -t14 * t32 + t8 * t17;
	t1 = t6 * t17 + t20 * t26;
	t2 = [0, t23 * t8 - t34 * t7, t24 * (t14 * t33 + t8 * t20) - t28 * t3, t3, (-t7 * t16 + t3 * t19) * r_i_i_C(1) + (-t3 * t16 - t7 * t19) * r_i_i_C(2), 0; 0, t23 * t6 - t34 * t5, t24 * (-t17 * t26 + t6 * t20) - t28 * t1, t1, (t1 * t19 - t5 * t16) * r_i_i_C(1) + (-t1 * t16 - t5 * t19) * r_i_i_C(2), 0; 1, (t23 * t18 + t34 * t21) * t15, -t28 * t9 + t24 * (t30 * t17 + t18 * t32), t9, (t16 * t31 + t9 * t19) * r_i_i_C(1) + (-t9 * t16 + t19 * t31) * r_i_i_C(2), 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:23:43
	% EndTime: 2019-10-09 22:23:43
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (197->46), mult. (509->82), div. (0->0), fcn. (652->10), ass. (0->36)
	t31 = sin(qJ(3));
	t34 = cos(qJ(3));
	t42 = pkin(3) + pkin(9) + r_i_i_C(2);
	t54 = qJ(4) * t31 + t42 * t34 + pkin(2);
	t53 = pkin(4) + pkin(8);
	t52 = pkin(5) + r_i_i_C(1);
	t29 = sin(pkin(6));
	t51 = t29 * t31;
	t50 = t29 * t34;
	t30 = sin(qJ(5));
	t49 = t30 * t31;
	t35 = cos(qJ(2));
	t48 = t30 * t35;
	t33 = cos(qJ(5));
	t47 = t31 * t33;
	t46 = t33 * t35;
	t45 = r_i_i_C(3) + qJ(6);
	t44 = cos(pkin(6));
	t43 = cos(pkin(10));
	t28 = sin(pkin(10));
	t40 = t28 * t44;
	t39 = t29 * t43;
	t38 = t44 * t43;
	t36 = t52 * t30 - t45 * t33 + qJ(4);
	t32 = sin(qJ(2));
	t23 = t32 * t51 - t44 * t34;
	t22 = -t32 * t40 + t43 * t35;
	t21 = t43 * t32 + t35 * t40;
	t20 = t28 * t35 + t32 * t38;
	t19 = t28 * t32 - t35 * t38;
	t15 = t23 * t33 + t29 * t48;
	t13 = t22 * t31 - t28 * t50;
	t11 = t20 * t31 + t34 * t39;
	t3 = -t13 * t33 + t21 * t30;
	t1 = -t11 * t33 + t19 * t30;
	t2 = [0, t52 * (-t21 * t49 + t22 * t33) + t45 * (t21 * t47 + t22 * t30) + t53 * t22 - t54 * t21, -t42 * t13 + t36 * (t22 * t34 + t28 * t51), t13, t45 * (t13 * t30 + t21 * t33) - t52 * t3, t3; 0, t52 * (-t19 * t49 + t20 * t33) + t45 * (t19 * t47 + t20 * t30) + t53 * t20 - t54 * t19, -t42 * t11 + t36 * (t20 * t34 - t31 * t39), t11, t45 * (t11 * t30 + t19 * t33) - t52 * t1, t1; 1, (t52 * (t31 * t48 + t32 * t33) + t45 * (t30 * t32 - t31 * t46) + t53 * t32 + t54 * t35) * t29, -t42 * t23 + t36 * (t44 * t31 + t32 * t50), t23, -t45 * (-t23 * t30 + t29 * t46) + t52 * t15, -t15;];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,6);
end