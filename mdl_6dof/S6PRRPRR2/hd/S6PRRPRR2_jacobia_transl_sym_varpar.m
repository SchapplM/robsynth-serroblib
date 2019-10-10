% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:27
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:28
	% EndTime: 2019-10-09 22:27:28
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:28
	% EndTime: 2019-10-09 22:27:28
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:28
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-09 22:27:28
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.15s
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
	% StartTime: 2019-10-09 22:27:28
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (70->28), mult. (123->52), div. (0->0), fcn. (150->10), ass. (0->24)
	t26 = r_i_i_C(3) + qJ(4) + pkin(8);
	t10 = sin(pkin(11));
	t11 = sin(pkin(6));
	t25 = t10 * t11;
	t12 = cos(pkin(11));
	t24 = t11 * t12;
	t16 = sin(qJ(2));
	t23 = t11 * t16;
	t17 = cos(qJ(3));
	t22 = t11 * t17;
	t13 = cos(pkin(6));
	t21 = t13 * t16;
	t18 = cos(qJ(2));
	t20 = t13 * t18;
	t9 = qJ(3) + pkin(12);
	t7 = sin(t9);
	t8 = cos(t9);
	t19 = t17 * pkin(3) + t8 * r_i_i_C(1) - t7 * r_i_i_C(2) + pkin(2);
	t15 = sin(qJ(3));
	t4 = -t10 * t21 + t12 * t18;
	t3 = t10 * t20 + t12 * t16;
	t2 = t10 * t18 + t12 * t21;
	t1 = t10 * t16 - t12 * t20;
	t5 = [0, -t19 * t3 + t26 * t4, (t8 * t25 - t4 * t7) * r_i_i_C(1) + (-t7 * t25 - t4 * t8) * r_i_i_C(2) + (t10 * t22 - t15 * t4) * pkin(3), t3, 0, 0; 0, -t19 * t1 + t26 * t2, (-t2 * t7 - t8 * t24) * r_i_i_C(1) + (-t2 * t8 + t7 * t24) * r_i_i_C(2) + (-t12 * t22 - t15 * t2) * pkin(3), t1, 0, 0; 1, (t26 * t16 + t19 * t18) * t11, (t13 * t8 - t7 * t23) * r_i_i_C(1) + (-t13 * t7 - t8 * t23) * r_i_i_C(2) + (t13 * t17 - t15 * t23) * pkin(3), -t11 * t18, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:28
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (176->42), mult. (300->74), div. (0->0), fcn. (376->12), ass. (0->32)
	t15 = qJ(3) + pkin(12);
	t13 = sin(t15);
	t14 = cos(t15);
	t25 = cos(qJ(3));
	t21 = sin(qJ(5));
	t24 = cos(qJ(5));
	t29 = t24 * r_i_i_C(1) - t21 * r_i_i_C(2) + pkin(4);
	t37 = pkin(9) + r_i_i_C(3);
	t38 = t25 * pkin(3) + t37 * t13 + t29 * t14 + pkin(2);
	t16 = sin(pkin(11));
	t17 = sin(pkin(6));
	t36 = t16 * t17;
	t18 = cos(pkin(11));
	t35 = t17 * t18;
	t23 = sin(qJ(2));
	t34 = t17 * t23;
	t33 = t17 * t25;
	t26 = cos(qJ(2));
	t32 = t17 * t26;
	t19 = cos(pkin(6));
	t31 = t19 * t23;
	t30 = t19 * t26;
	t28 = t21 * r_i_i_C(1) + t24 * r_i_i_C(2) + pkin(8) + qJ(4);
	t22 = sin(qJ(3));
	t10 = -t16 * t31 + t18 * t26;
	t9 = t16 * t30 + t18 * t23;
	t8 = t16 * t26 + t18 * t31;
	t7 = t16 * t23 - t18 * t30;
	t6 = t19 * t13 + t14 * t34;
	t4 = t10 * t14 + t13 * t36;
	t2 = -t13 * t35 + t8 * t14;
	t1 = [0, t28 * t10 - t38 * t9, t37 * t4 + (-t10 * t22 + t16 * t33) * pkin(3) + t29 * (-t10 * t13 + t14 * t36), t9, (-t4 * t21 + t9 * t24) * r_i_i_C(1) + (-t9 * t21 - t4 * t24) * r_i_i_C(2), 0; 0, t28 * t8 - t38 * t7, t37 * t2 + (-t18 * t33 - t22 * t8) * pkin(3) + t29 * (-t8 * t13 - t14 * t35), t7, (-t2 * t21 + t7 * t24) * r_i_i_C(1) + (-t2 * t24 - t7 * t21) * r_i_i_C(2), 0; 1, (t28 * t23 + t38 * t26) * t17, t37 * t6 + (t19 * t25 - t22 * t34) * pkin(3) + t29 * (-t13 * t34 + t19 * t14), -t32, (-t6 * t21 - t24 * t32) * r_i_i_C(1) + (t21 * t32 - t6 * t24) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:28
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (289->51), mult. (415->85), div. (0->0), fcn. (520->14), ass. (0->38)
	t24 = qJ(3) + pkin(12);
	t20 = sin(t24);
	t21 = cos(t24);
	t35 = cos(qJ(3));
	t25 = qJ(5) + qJ(6);
	t22 = sin(t25);
	t23 = cos(t25);
	t34 = cos(qJ(5));
	t40 = t34 * pkin(5) + r_i_i_C(1) * t23 - r_i_i_C(2) * t22 + pkin(4);
	t48 = r_i_i_C(3) + pkin(10) + pkin(9);
	t52 = t35 * pkin(3) + t48 * t20 + t40 * t21 + pkin(2);
	t26 = sin(pkin(11));
	t28 = cos(pkin(11));
	t33 = sin(qJ(2));
	t29 = cos(pkin(6));
	t36 = cos(qJ(2));
	t41 = t29 * t36;
	t13 = t26 * t33 - t28 * t41;
	t42 = t29 * t33;
	t14 = t26 * t36 + t28 * t42;
	t27 = sin(pkin(6));
	t46 = t27 * t28;
	t8 = t14 * t21 - t20 * t46;
	t51 = (t13 * t23 - t8 * t22) * r_i_i_C(1) + (-t13 * t22 - t8 * t23) * r_i_i_C(2);
	t16 = -t26 * t42 + t28 * t36;
	t47 = t26 * t27;
	t10 = t16 * t21 + t20 * t47;
	t15 = t26 * t41 + t28 * t33;
	t50 = (-t10 * t22 + t15 * t23) * r_i_i_C(1) + (-t10 * t23 - t15 * t22) * r_i_i_C(2);
	t45 = t27 * t33;
	t12 = t29 * t20 + t21 * t45;
	t43 = t27 * t36;
	t49 = (-t12 * t22 - t23 * t43) * r_i_i_C(1) + (-t12 * t23 + t22 * t43) * r_i_i_C(2);
	t44 = t27 * t35;
	t31 = sin(qJ(5));
	t39 = t31 * pkin(5) + t22 * r_i_i_C(1) + t23 * r_i_i_C(2) + pkin(8) + qJ(4);
	t32 = sin(qJ(3));
	t1 = [0, -t15 * t52 + t16 * t39, t48 * t10 + (-t16 * t32 + t26 * t44) * pkin(3) + t40 * (-t16 * t20 + t21 * t47), t15, (-t10 * t31 + t15 * t34) * pkin(5) + t50, t50; 0, -t13 * t52 + t14 * t39, t48 * t8 + (-t14 * t32 - t28 * t44) * pkin(3) + t40 * (-t14 * t20 - t21 * t46), t13, (t13 * t34 - t31 * t8) * pkin(5) + t51, t51; 1, (t39 * t33 + t52 * t36) * t27, t48 * t12 + (t29 * t35 - t32 * t45) * pkin(3) + t40 * (-t20 * t45 + t29 * t21), -t43, (-t12 * t31 - t34 * t43) * pkin(5) + t49, t49;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end