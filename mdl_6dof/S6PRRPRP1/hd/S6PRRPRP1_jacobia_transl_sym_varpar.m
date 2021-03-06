% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRP1
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
% Datum: 2019-10-09 22:16
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRP1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRP1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:16:24
	% EndTime: 2019-10-09 22:16:24
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:16:24
	% EndTime: 2019-10-09 22:16:24
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:16:24
	% EndTime: 2019-10-09 22:16:24
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-10-09 22:16:24
	% EndTime: 2019-10-09 22:16:24
	% DurationCPUTime: 0.15s
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
	% StartTime: 2019-10-09 22:16:24
	% EndTime: 2019-10-09 22:16:24
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (70->28), mult. (123->52), div. (0->0), fcn. (150->10), ass. (0->24)
	t26 = r_i_i_C(3) + qJ(4) + pkin(8);
	t10 = sin(pkin(10));
	t11 = sin(pkin(6));
	t25 = t10 * t11;
	t12 = cos(pkin(10));
	t24 = t11 * t12;
	t16 = sin(qJ(2));
	t23 = t11 * t16;
	t17 = cos(qJ(3));
	t22 = t11 * t17;
	t13 = cos(pkin(6));
	t21 = t13 * t16;
	t18 = cos(qJ(2));
	t20 = t13 * t18;
	t9 = qJ(3) + pkin(11);
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
	% StartTime: 2019-10-09 22:16:24
	% EndTime: 2019-10-09 22:16:24
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (176->42), mult. (300->74), div. (0->0), fcn. (376->12), ass. (0->32)
	t15 = qJ(3) + pkin(11);
	t13 = sin(t15);
	t14 = cos(t15);
	t25 = cos(qJ(3));
	t21 = sin(qJ(5));
	t24 = cos(qJ(5));
	t29 = t24 * r_i_i_C(1) - t21 * r_i_i_C(2) + pkin(4);
	t37 = pkin(9) + r_i_i_C(3);
	t38 = t25 * pkin(3) + t37 * t13 + t29 * t14 + pkin(2);
	t16 = sin(pkin(10));
	t17 = sin(pkin(6));
	t36 = t16 * t17;
	t18 = cos(pkin(10));
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
	% StartTime: 2019-10-09 22:16:24
	% EndTime: 2019-10-09 22:16:24
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (223->43), mult. (368->74), div. (0->0), fcn. (461->12), ass. (0->36)
	t46 = pkin(5) + r_i_i_C(1);
	t18 = qJ(3) + pkin(11);
	t16 = sin(t18);
	t17 = cos(t18);
	t28 = cos(qJ(3));
	t24 = sin(qJ(5));
	t27 = cos(qJ(5));
	t33 = -t24 * r_i_i_C(2) + t46 * t27 + pkin(4);
	t44 = r_i_i_C(3) + qJ(6) + pkin(9);
	t45 = t28 * pkin(3) + t44 * t16 + t33 * t17 + pkin(2);
	t19 = sin(pkin(10));
	t20 = sin(pkin(6));
	t43 = t19 * t20;
	t21 = cos(pkin(10));
	t42 = t20 * t21;
	t26 = sin(qJ(2));
	t41 = t20 * t26;
	t40 = t20 * t28;
	t29 = cos(qJ(2));
	t39 = t20 * t29;
	t38 = cos(pkin(6));
	t37 = t26 * t38;
	t36 = t29 * t38;
	t31 = t27 * r_i_i_C(2) + t46 * t24 + pkin(8) + qJ(4);
	t25 = sin(qJ(3));
	t10 = -t19 * t37 + t21 * t29;
	t9 = t19 * t36 + t21 * t26;
	t8 = t19 * t29 + t21 * t37;
	t7 = t19 * t26 - t21 * t36;
	t6 = t38 * t16 + t17 * t41;
	t5 = t16 * t41 - t38 * t17;
	t4 = t10 * t17 + t16 * t43;
	t3 = t10 * t16 - t17 * t43;
	t2 = -t16 * t42 + t8 * t17;
	t1 = t8 * t16 + t17 * t42;
	t11 = [0, t31 * t10 - t45 * t9, t44 * t4 + (-t10 * t25 + t19 * t40) * pkin(3) - t33 * t3, t9, (-t9 * t24 - t4 * t27) * r_i_i_C(2) + t46 * (-t4 * t24 + t9 * t27), t3; 0, t31 * t8 - t45 * t7, t44 * t2 + (-t21 * t40 - t25 * t8) * pkin(3) - t33 * t1, t7, (-t2 * t27 - t7 * t24) * r_i_i_C(2) + t46 * (-t2 * t24 + t7 * t27), t1; 1, (t31 * t26 + t29 * t45) * t20, t44 * t6 + (-t25 * t41 + t38 * t28) * pkin(3) - t33 * t5, -t39, (t24 * t39 - t6 * t27) * r_i_i_C(2) + t46 * (-t6 * t24 - t27 * t39), t5;];
	Ja_transl = t11;
else
	Ja_transl=NaN(3,6);
end