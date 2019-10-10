% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:20
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRP3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRP3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:20:05
	% EndTime: 2019-10-09 22:20:05
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:20:05
	% EndTime: 2019-10-09 22:20:05
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:20:05
	% EndTime: 2019-10-09 22:20:05
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-09 22:20:05
	% EndTime: 2019-10-09 22:20:05
	% DurationCPUTime: 0.13s
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
	% StartTime: 2019-10-09 22:20:05
	% EndTime: 2019-10-09 22:20:05
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (84->23), mult. (222->43), div. (0->0), fcn. (279->10), ass. (0->24)
	t15 = sin(pkin(6));
	t19 = sin(qJ(3));
	t30 = t15 * t19;
	t21 = cos(qJ(3));
	t29 = t15 * t21;
	t18 = cos(pkin(6));
	t20 = sin(qJ(2));
	t28 = t18 * t20;
	t22 = cos(qJ(2));
	t27 = t18 * t22;
	t26 = r_i_i_C(3) + qJ(4);
	t13 = sin(pkin(11));
	t16 = cos(pkin(11));
	t25 = r_i_i_C(1) * t16 - r_i_i_C(2) * t13 + pkin(3);
	t24 = t13 * r_i_i_C(1) + t16 * r_i_i_C(2) + pkin(8);
	t23 = t26 * t19 + t25 * t21 + pkin(2);
	t17 = cos(pkin(10));
	t14 = sin(pkin(10));
	t9 = -t18 * t21 + t20 * t30;
	t8 = -t14 * t28 + t17 * t22;
	t6 = t14 * t22 + t17 * t28;
	t3 = -t14 * t29 + t19 * t8;
	t1 = t17 * t29 + t19 * t6;
	t2 = [0, t24 * t8 + t23 * (-t14 * t27 - t17 * t20), t26 * (t14 * t30 + t21 * t8) - t25 * t3, t3, 0, 0; 0, t24 * t6 + t23 * (-t14 * t20 + t17 * t27), t26 * (-t17 * t30 + t21 * t6) - t25 * t1, t1, 0, 0; 1, (t24 * t20 + t23 * t22) * t15, t26 * (t18 * t19 + t20 * t29) - t25 * t9, t9, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:20:05
	% EndTime: 2019-10-09 22:20:05
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (154->35), mult. (303->65), div. (0->0), fcn. (382->12), ass. (0->32)
	t22 = sin(qJ(3));
	t24 = cos(qJ(3));
	t17 = pkin(11) + qJ(5);
	t15 = sin(t17);
	t16 = cos(t17);
	t28 = t16 * r_i_i_C(1) - t15 * r_i_i_C(2) + cos(pkin(11)) * pkin(4) + pkin(3);
	t37 = r_i_i_C(3) + pkin(9) + qJ(4);
	t38 = t37 * t22 + t28 * t24 + pkin(2);
	t20 = sin(pkin(6));
	t36 = t20 * t22;
	t35 = t20 * t24;
	t25 = cos(qJ(2));
	t34 = t20 * t25;
	t33 = cos(pkin(6));
	t32 = cos(pkin(10));
	t19 = sin(pkin(10));
	t31 = t19 * t33;
	t30 = t20 * t32;
	t29 = t33 * t32;
	t27 = sin(pkin(11)) * pkin(4) + t15 * r_i_i_C(1) + t16 * r_i_i_C(2) + pkin(8);
	t23 = sin(qJ(2));
	t10 = t33 * t22 + t23 * t35;
	t9 = t23 * t36 - t33 * t24;
	t8 = -t23 * t31 + t32 * t25;
	t7 = t32 * t23 + t25 * t31;
	t6 = t19 * t25 + t23 * t29;
	t5 = t19 * t23 - t25 * t29;
	t4 = t19 * t36 + t8 * t24;
	t3 = -t19 * t35 + t8 * t22;
	t2 = -t22 * t30 + t6 * t24;
	t1 = t6 * t22 + t24 * t30;
	t11 = [0, t27 * t8 - t38 * t7, -t28 * t3 + t37 * t4, t3, (-t4 * t15 + t7 * t16) * r_i_i_C(1) + (-t7 * t15 - t4 * t16) * r_i_i_C(2), 0; 0, t27 * t6 - t38 * t5, -t28 * t1 + t37 * t2, t1, (-t2 * t15 + t5 * t16) * r_i_i_C(1) + (-t5 * t15 - t2 * t16) * r_i_i_C(2), 0; 1, (t27 * t23 + t38 * t25) * t20, t37 * t10 - t28 * t9, t9, (-t10 * t15 - t16 * t34) * r_i_i_C(1) + (-t10 * t16 + t15 * t34) * r_i_i_C(2), 0;];
	Ja_transl = t11;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:20:05
	% EndTime: 2019-10-09 22:20:05
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (262->47), mult. (487->84), div. (0->0), fcn. (623->12), ass. (0->41)
	t28 = cos(pkin(11)) * pkin(4) + pkin(3);
	t36 = sin(qJ(3));
	t38 = cos(qJ(3));
	t56 = r_i_i_C(2) + pkin(9) + qJ(4);
	t58 = t28 * t38 + t56 * t36 + pkin(2);
	t57 = pkin(5) + r_i_i_C(1);
	t31 = pkin(11) + qJ(5);
	t29 = sin(t31);
	t55 = t29 * t38;
	t30 = cos(t31);
	t54 = t30 * t38;
	t34 = sin(pkin(6));
	t53 = t34 * t36;
	t52 = t34 * t38;
	t39 = cos(qJ(2));
	t51 = t34 * t39;
	t50 = t38 * t39;
	t49 = r_i_i_C(3) + qJ(6);
	t48 = cos(pkin(6));
	t47 = cos(pkin(10));
	t46 = pkin(4) * sin(pkin(11)) + pkin(8);
	t33 = sin(pkin(10));
	t44 = t33 * t48;
	t43 = t34 * t47;
	t42 = t48 * t47;
	t40 = -t49 * t29 - t57 * t30 - t28;
	t37 = sin(qJ(2));
	t24 = t48 * t36 + t37 * t52;
	t23 = t37 * t53 - t48 * t38;
	t22 = -t37 * t44 + t47 * t39;
	t21 = t47 * t37 + t39 * t44;
	t20 = t33 * t39 + t37 * t42;
	t19 = t33 * t37 - t39 * t42;
	t14 = t22 * t38 + t33 * t53;
	t13 = t22 * t36 - t33 * t52;
	t12 = t20 * t38 - t36 * t43;
	t11 = t20 * t36 + t38 * t43;
	t9 = t24 * t29 + t30 * t51;
	t3 = t14 * t29 - t21 * t30;
	t1 = t12 * t29 - t19 * t30;
	t2 = [0, t57 * (-t21 * t54 + t22 * t29) + t49 * (-t21 * t55 - t22 * t30) + t46 * t22 - t58 * t21, t40 * t13 + t56 * t14, t13, t49 * (t14 * t30 + t21 * t29) - t57 * t3, t3; 0, t57 * (-t19 * t54 + t20 * t29) + t49 * (-t19 * t55 - t20 * t30) + t46 * t20 - t58 * t19, t40 * t11 + t56 * t12, t11, t49 * (t12 * t30 + t19 * t29) - t57 * t1, t1; 1, (t57 * (t29 * t37 + t30 * t50) + t49 * (t29 * t50 - t30 * t37) + t46 * t37 + t58 * t39) * t34, t40 * t23 + t56 * t24, t23, -t57 * t9 + t49 * (t24 * t30 - t29 * t51), t9;];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,6);
end