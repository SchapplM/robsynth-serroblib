% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
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
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (25->11), mult. (63->20), div. (0->0), fcn. (78->8), ass. (0->13)
	t11 = cos(pkin(6));
	t12 = sin(qJ(2));
	t17 = t11 * t12;
	t13 = cos(qJ(2));
	t16 = t11 * t13;
	t15 = r_i_i_C(3) + qJ(3);
	t14 = r_i_i_C(1) * cos(pkin(11)) - r_i_i_C(2) * sin(pkin(11)) + pkin(2);
	t10 = cos(pkin(10));
	t8 = sin(pkin(6));
	t7 = sin(pkin(10));
	t3 = t10 * t12 + t7 * t16;
	t1 = -t10 * t16 + t7 * t12;
	t2 = [0, t15 * (t10 * t13 - t7 * t17) - t14 * t3, t3, 0, 0, 0; 0, t15 * (t10 * t17 + t7 * t13) - t14 * t1, t1, 0, 0, 0; 1, (t15 * t12 + t14 * t13) * t8, -t8 * t13, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (62->22), mult. (102->42), div. (0->0), fcn. (127->9), ass. (0->21)
	t23 = r_i_i_C(3) + pkin(8) + qJ(3);
	t10 = sin(pkin(10));
	t11 = sin(pkin(6));
	t22 = t10 * t11;
	t12 = cos(pkin(10));
	t21 = t11 * t12;
	t15 = sin(qJ(2));
	t20 = t11 * t15;
	t13 = cos(pkin(6));
	t19 = t13 * t15;
	t16 = cos(qJ(2));
	t18 = t13 * t16;
	t9 = pkin(11) + qJ(4);
	t7 = sin(t9);
	t8 = cos(t9);
	t17 = r_i_i_C(1) * t8 - r_i_i_C(2) * t7 + cos(pkin(11)) * pkin(3) + pkin(2);
	t4 = -t10 * t19 + t12 * t16;
	t3 = t10 * t18 + t12 * t15;
	t2 = t10 * t16 + t12 * t19;
	t1 = t10 * t15 - t12 * t18;
	t5 = [0, -t17 * t3 + t23 * t4, t3, (t8 * t22 - t4 * t7) * r_i_i_C(1) + (-t7 * t22 - t4 * t8) * r_i_i_C(2), 0, 0; 0, -t17 * t1 + t23 * t2, t1, (-t2 * t7 - t8 * t21) * r_i_i_C(1) + (-t2 * t8 + t7 * t21) * r_i_i_C(2), 0, 0; 1, (t23 * t15 + t17 * t16) * t11, -t11 * t16, (t13 * t8 - t7 * t20) * r_i_i_C(1) + (-t13 * t7 - t8 * t20) * r_i_i_C(2), 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (117->24), mult. (179->42), div. (0->0), fcn. (226->9), ass. (0->26)
	t15 = pkin(11) + qJ(4);
	t13 = sin(t15);
	t14 = cos(t15);
	t24 = r_i_i_C(3) + qJ(5);
	t31 = pkin(4) - r_i_i_C(2);
	t32 = t24 * t13 + t31 * t14 + cos(pkin(11)) * pkin(3) + pkin(2);
	t30 = r_i_i_C(1) + pkin(8) + qJ(3);
	t16 = sin(pkin(10));
	t17 = sin(pkin(6));
	t29 = t16 * t17;
	t18 = cos(pkin(10));
	t28 = t17 * t18;
	t21 = sin(qJ(2));
	t27 = t17 * t21;
	t19 = cos(pkin(6));
	t26 = t19 * t21;
	t22 = cos(qJ(2));
	t25 = t19 * t22;
	t10 = -t16 * t26 + t18 * t22;
	t9 = t16 * t25 + t18 * t21;
	t8 = t16 * t22 + t18 * t26;
	t7 = t16 * t21 - t18 * t25;
	t5 = t13 * t27 - t14 * t19;
	t3 = t10 * t13 - t14 * t29;
	t1 = t13 * t8 + t14 * t28;
	t2 = [0, t30 * t10 - t32 * t9, t9, t24 * (t10 * t14 + t13 * t29) - t31 * t3, t3, 0; 0, t30 * t8 - t32 * t7, t7, t24 * (-t13 * t28 + t14 * t8) - t31 * t1, t1, 0; 1, (t30 * t21 + t32 * t22) * t17, -t17 * t22, t24 * (t13 * t19 + t14 * t27) - t31 * t5, t5, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (206->37), mult. (335->64), div. (0->0), fcn. (426->11), ass. (0->29)
	t17 = pkin(11) + qJ(4);
	t15 = sin(t17);
	t16 = cos(t17);
	t21 = sin(qJ(6));
	t23 = cos(qJ(6));
	t27 = t21 * r_i_i_C(1) + t23 * r_i_i_C(2) + qJ(5);
	t31 = pkin(4) + pkin(9) + r_i_i_C(3);
	t37 = t27 * t15 + t31 * t16 + cos(pkin(11)) * pkin(3) + pkin(2);
	t18 = sin(pkin(10));
	t19 = sin(pkin(6));
	t36 = t18 * t19;
	t22 = sin(qJ(2));
	t35 = t19 * t22;
	t24 = cos(qJ(2));
	t34 = t19 * t24;
	t33 = cos(pkin(6));
	t32 = cos(pkin(10));
	t30 = t18 * t33;
	t29 = t19 * t32;
	t28 = t33 * t32;
	t26 = t23 * r_i_i_C(1) - t21 * r_i_i_C(2) + pkin(5) + pkin(8) + qJ(3);
	t10 = -t22 * t30 + t32 * t24;
	t9 = t32 * t22 + t24 * t30;
	t8 = t18 * t24 + t22 * t28;
	t7 = t18 * t22 - t24 * t28;
	t5 = t15 * t35 - t33 * t16;
	t3 = t10 * t15 - t16 * t36;
	t1 = t8 * t15 + t16 * t29;
	t2 = [0, t26 * t10 - t37 * t9, t9, t27 * (t10 * t16 + t15 * t36) - t31 * t3, t3, (-t9 * t21 + t3 * t23) * r_i_i_C(1) + (-t3 * t21 - t9 * t23) * r_i_i_C(2); 0, t26 * t8 - t37 * t7, t7, t27 * (-t15 * t29 + t8 * t16) - t31 * t1, t1, (t1 * t23 - t7 * t21) * r_i_i_C(1) + (-t1 * t21 - t7 * t23) * r_i_i_C(2); 1, (t26 * t22 + t37 * t24) * t19, -t34, t27 * (t33 * t15 + t16 * t35) - t31 * t5, t5, (t21 * t34 + t5 * t23) * r_i_i_C(1) + (-t5 * t21 + t23 * t34) * r_i_i_C(2);];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,6);
end