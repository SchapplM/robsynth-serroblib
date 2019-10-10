% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPPR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.12s
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
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.12s
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
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (66->21), mult. (169->39), div. (0->0), fcn. (211->8), ass. (0->22)
	t28 = pkin(3) + r_i_i_C(1);
	t27 = pkin(8) + r_i_i_C(2);
	t14 = sin(pkin(6));
	t17 = sin(qJ(3));
	t26 = t14 * t17;
	t19 = cos(qJ(3));
	t25 = t14 * t19;
	t16 = cos(pkin(6));
	t18 = sin(qJ(2));
	t24 = t16 * t18;
	t20 = cos(qJ(2));
	t23 = t16 * t20;
	t22 = r_i_i_C(3) + qJ(4);
	t21 = t22 * t17 + t19 * t28 + pkin(2);
	t15 = cos(pkin(10));
	t13 = sin(pkin(10));
	t9 = -t16 * t19 + t18 * t26;
	t8 = -t13 * t24 + t15 * t20;
	t6 = t13 * t20 + t15 * t24;
	t3 = -t13 * t25 + t17 * t8;
	t1 = t15 * t25 + t17 * t6;
	t2 = [0, t27 * t8 + t21 * (-t13 * t23 - t15 * t18), t22 * (t13 * t26 + t19 * t8) - t28 * t3, t3, 0, 0; 0, t27 * t6 + t21 * (-t13 * t18 + t15 * t23), t22 * (-t15 * t26 + t19 * t6) - t28 * t1, t1, 0, 0; 1, (t18 * t27 + t20 * t21) * t14, -t28 * t9 + t22 * (t16 * t17 + t18 * t25), t9, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (86->21), mult. (217->40), div. (0->0), fcn. (273->8), ass. (0->24)
	t15 = sin(pkin(6));
	t18 = sin(qJ(3));
	t29 = t15 * t18;
	t20 = cos(qJ(3));
	t28 = t15 * t20;
	t17 = cos(pkin(6));
	t19 = sin(qJ(2));
	t27 = t17 * t19;
	t21 = cos(qJ(2));
	t26 = t17 * t21;
	t25 = r_i_i_C(1) + qJ(4);
	t24 = pkin(3) + pkin(4) - r_i_i_C(2);
	t23 = pkin(8) - r_i_i_C(3) - qJ(5);
	t22 = t25 * t18 + t24 * t20 + pkin(2);
	t16 = cos(pkin(10));
	t14 = sin(pkin(10));
	t9 = -t17 * t20 + t19 * t29;
	t8 = -t14 * t27 + t16 * t21;
	t7 = -t14 * t26 - t16 * t19;
	t6 = t14 * t21 + t16 * t27;
	t5 = -t14 * t19 + t16 * t26;
	t3 = -t14 * t28 + t18 * t8;
	t1 = t16 * t28 + t18 * t6;
	t2 = [0, t22 * t7 + t23 * t8, t25 * (t14 * t29 + t20 * t8) - t24 * t3, t3, t7, 0; 0, t22 * t5 + t23 * t6, t25 * (-t16 * t29 + t20 * t6) - t24 * t1, t1, t5, 0; 1, (t23 * t19 + t22 * t21) * t15, t25 * (t17 * t18 + t19 * t28) - t24 * t9, t9, t15 * t21, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (153->34), mult. (394->62), div. (0->0), fcn. (499->10), ass. (0->27)
	t15 = sin(pkin(6));
	t19 = sin(qJ(3));
	t32 = t15 * t19;
	t22 = cos(qJ(3));
	t31 = t15 * t22;
	t23 = cos(qJ(2));
	t30 = t15 * t23;
	t17 = cos(pkin(6));
	t20 = sin(qJ(2));
	t29 = t17 * t20;
	t28 = t17 * t23;
	t27 = pkin(3) + pkin(4) + pkin(9) + r_i_i_C(3);
	t18 = sin(qJ(6));
	t21 = cos(qJ(6));
	t26 = t21 * r_i_i_C(1) - t18 * r_i_i_C(2) + pkin(5) + qJ(4);
	t25 = -t18 * r_i_i_C(1) - t21 * r_i_i_C(2) + pkin(8) - qJ(5);
	t24 = t26 * t19 + t27 * t22 + pkin(2);
	t16 = cos(pkin(10));
	t14 = sin(pkin(10));
	t9 = -t17 * t22 + t20 * t32;
	t8 = -t14 * t29 + t16 * t23;
	t7 = -t14 * t28 - t16 * t20;
	t6 = t14 * t23 + t16 * t29;
	t5 = -t14 * t20 + t16 * t28;
	t3 = -t14 * t31 + t8 * t19;
	t1 = t16 * t31 + t6 * t19;
	t2 = [0, t24 * t7 + t25 * t8, t26 * (t14 * t32 + t8 * t22) - t27 * t3, t3, t7, (-t3 * t18 + t7 * t21) * r_i_i_C(1) + (-t7 * t18 - t3 * t21) * r_i_i_C(2); 0, t24 * t5 + t25 * t6, t26 * (-t16 * t32 + t6 * t22) - t27 * t1, t1, t5, (-t1 * t18 + t5 * t21) * r_i_i_C(1) + (-t1 * t21 - t5 * t18) * r_i_i_C(2); 1, (t25 * t20 + t24 * t23) * t15, -t27 * t9 + t26 * (t17 * t19 + t20 * t31), t9, t30, (-t9 * t18 + t21 * t30) * r_i_i_C(1) + (-t18 * t30 - t9 * t21) * r_i_i_C(2);];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,6);
end