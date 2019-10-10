% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:14
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPPR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
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
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
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
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
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
	t19 = t20 * t15 + t17 * t26 + pkin(2);
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
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (107->24), mult. (278->43), div. (0->0), fcn. (352->10), ass. (0->27)
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
	t26 = pkin(3) + r_i_i_C(3) + qJ(5);
	t13 = sin(pkin(11));
	t16 = cos(pkin(11));
	t25 = r_i_i_C(1) * t13 + r_i_i_C(2) * t16 + qJ(4);
	t24 = t16 * r_i_i_C(1) - t13 * r_i_i_C(2) + pkin(4) + pkin(8);
	t23 = t25 * t19 + t26 * t21 + pkin(2);
	t17 = cos(pkin(10));
	t14 = sin(pkin(10));
	t10 = t18 * t19 + t20 * t29;
	t9 = -t18 * t21 + t20 * t30;
	t8 = -t14 * t28 + t17 * t22;
	t6 = t14 * t22 + t17 * t28;
	t4 = t14 * t30 + t8 * t21;
	t3 = -t14 * t29 + t8 * t19;
	t2 = -t17 * t30 + t6 * t21;
	t1 = t17 * t29 + t6 * t19;
	t5 = [0, t24 * t8 + t23 * (-t14 * t27 - t17 * t20), t25 * t4 - t26 * t3, t3, t4, 0; 0, t24 * t6 + t23 * (-t14 * t20 + t17 * t27), -t26 * t1 + t25 * t2, t1, t2, 0; 1, (t24 * t20 + t23 * t22) * t15, t25 * t10 - t26 * t9, t9, t10, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (182->36), mult. (380->65), div. (0->0), fcn. (481->12), ass. (0->32)
	t22 = sin(qJ(3));
	t24 = cos(qJ(3));
	t17 = pkin(11) + qJ(6);
	t15 = sin(t17);
	t16 = cos(t17);
	t27 = pkin(5) * sin(pkin(11)) + t15 * r_i_i_C(1) + t16 * r_i_i_C(2) + qJ(4);
	t32 = pkin(3) + r_i_i_C(3) + pkin(9) + qJ(5);
	t38 = t27 * t22 + t32 * t24 + pkin(2);
	t20 = sin(pkin(6));
	t37 = t20 * t22;
	t36 = t20 * t24;
	t25 = cos(qJ(2));
	t35 = t20 * t25;
	t34 = cos(pkin(6));
	t33 = cos(pkin(10));
	t19 = sin(pkin(10));
	t31 = t19 * t34;
	t30 = t20 * t33;
	t29 = t34 * t33;
	t28 = t16 * r_i_i_C(1) - t15 * r_i_i_C(2) + pkin(8) + cos(pkin(11)) * pkin(5) + pkin(4);
	t23 = sin(qJ(2));
	t10 = t34 * t22 + t23 * t36;
	t9 = t23 * t37 - t34 * t24;
	t8 = -t23 * t31 + t33 * t25;
	t7 = t33 * t23 + t25 * t31;
	t6 = t19 * t25 + t23 * t29;
	t5 = t19 * t23 - t25 * t29;
	t4 = t19 * t37 + t8 * t24;
	t3 = -t19 * t36 + t8 * t22;
	t2 = -t22 * t30 + t6 * t24;
	t1 = t6 * t22 + t24 * t30;
	t11 = [0, t28 * t8 - t38 * t7, t27 * t4 - t32 * t3, t3, t4, (-t7 * t15 + t3 * t16) * r_i_i_C(1) + (-t3 * t15 - t7 * t16) * r_i_i_C(2); 0, t28 * t6 - t38 * t5, -t32 * t1 + t27 * t2, t1, t2, (t1 * t16 - t5 * t15) * r_i_i_C(1) + (-t1 * t15 - t5 * t16) * r_i_i_C(2); 1, (t28 * t23 + t38 * t25) * t20, t27 * t10 - t32 * t9, t9, t10, (t15 * t35 + t9 * t16) * r_i_i_C(1) + (-t9 * t15 + t16 * t35) * r_i_i_C(2);];
	Ja_transl = t11;
else
	Ja_transl=NaN(3,6);
end