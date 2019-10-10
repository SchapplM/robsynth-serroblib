% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
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
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
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
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
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
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (144->27), mult. (232->46), div. (0->0), fcn. (294->11), ass. (0->28)
	t17 = pkin(11) + qJ(4);
	t15 = sin(t17);
	t16 = cos(t17);
	t18 = sin(pkin(12));
	t21 = cos(pkin(12));
	t27 = r_i_i_C(1) * t21 - r_i_i_C(2) * t18 + pkin(4);
	t33 = r_i_i_C(3) + qJ(5);
	t36 = t33 * t15 + t27 * t16 + cos(pkin(11)) * pkin(3) + pkin(2);
	t19 = sin(pkin(10));
	t20 = sin(pkin(6));
	t35 = t19 * t20;
	t23 = sin(qJ(2));
	t34 = t20 * t23;
	t32 = cos(pkin(6));
	t31 = cos(pkin(10));
	t30 = t19 * t32;
	t29 = t20 * t31;
	t28 = t32 * t31;
	t26 = t18 * r_i_i_C(1) + t21 * r_i_i_C(2) + pkin(8) + qJ(3);
	t24 = cos(qJ(2));
	t10 = -t23 * t30 + t31 * t24;
	t9 = t31 * t23 + t24 * t30;
	t8 = t19 * t24 + t23 * t28;
	t7 = t19 * t23 - t24 * t28;
	t5 = t15 * t34 - t32 * t16;
	t3 = t10 * t15 - t16 * t35;
	t1 = t8 * t15 + t16 * t29;
	t2 = [0, t26 * t10 - t36 * t9, t9, t33 * (t10 * t16 + t15 * t35) - t27 * t3, t3, 0; 0, t26 * t8 - t36 * t7, t7, t33 * (-t15 * t29 + t8 * t16) - t27 * t1, t1, 0; 1, (t26 * t23 + t36 * t24) * t20, -t20 * t24, t33 * (t32 * t15 + t16 * t34) - t27 * t5, t5, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (226->39), mult. (313->66), div. (0->0), fcn. (397->13), ass. (0->33)
	t21 = pkin(11) + qJ(4);
	t17 = sin(t21);
	t19 = cos(t21);
	t20 = pkin(12) + qJ(6);
	t16 = sin(t20);
	t18 = cos(t20);
	t31 = t18 * r_i_i_C(1) - t16 * r_i_i_C(2) + cos(pkin(12)) * pkin(5) + pkin(4);
	t40 = r_i_i_C(3) + pkin(9) + qJ(5);
	t41 = t40 * t17 + t31 * t19 + cos(pkin(11)) * pkin(3) + pkin(2);
	t23 = sin(pkin(10));
	t24 = sin(pkin(6));
	t39 = t23 * t24;
	t27 = sin(qJ(2));
	t38 = t24 * t27;
	t28 = cos(qJ(2));
	t37 = t24 * t28;
	t36 = cos(pkin(6));
	t35 = cos(pkin(10));
	t34 = t23 * t36;
	t33 = t24 * t35;
	t32 = t36 * t35;
	t30 = sin(pkin(12)) * pkin(5) + t16 * r_i_i_C(1) + t18 * r_i_i_C(2) + pkin(8) + qJ(3);
	t10 = -t27 * t34 + t35 * t28;
	t9 = t35 * t27 + t28 * t34;
	t8 = t23 * t28 + t27 * t32;
	t7 = t23 * t27 - t28 * t32;
	t6 = t36 * t17 + t19 * t38;
	t5 = t17 * t38 - t36 * t19;
	t4 = t10 * t19 + t17 * t39;
	t3 = t10 * t17 - t19 * t39;
	t2 = -t17 * t33 + t8 * t19;
	t1 = t8 * t17 + t19 * t33;
	t11 = [0, t30 * t10 - t41 * t9, t9, -t31 * t3 + t40 * t4, t3, (-t4 * t16 + t9 * t18) * r_i_i_C(1) + (-t9 * t16 - t4 * t18) * r_i_i_C(2); 0, t30 * t8 - t41 * t7, t7, -t31 * t1 + t40 * t2, t1, (-t2 * t16 + t7 * t18) * r_i_i_C(1) + (-t7 * t16 - t2 * t18) * r_i_i_C(2); 1, (t30 * t27 + t41 * t28) * t24, -t37, -t31 * t5 + t40 * t6, t5, (-t6 * t16 - t18 * t37) * r_i_i_C(1) + (t16 * t37 - t6 * t18) * r_i_i_C(2);];
	Ja_transl = t11;
else
	Ja_transl=NaN(3,6);
end