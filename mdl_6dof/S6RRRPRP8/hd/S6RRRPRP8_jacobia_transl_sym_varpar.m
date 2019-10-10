% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:47
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP8_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP8_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP8_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:47:48
	% EndTime: 2019-10-10 11:47:48
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:47:48
	% EndTime: 2019-10-10 11:47:48
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:47:48
	% EndTime: 2019-10-10 11:47:48
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (11->6), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->9)
	t8 = pkin(7) + r_i_i_C(3);
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = t3 * r_i_i_C(1) - t1 * r_i_i_C(2);
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t5 * t2 + t8 * t4, t6 * t4, 0, 0, 0, 0; t8 * t2 + t5 * t4, t6 * t2, 0, 0, 0, 0; 0, t7, 0, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:47:48
	% EndTime: 2019-10-10 11:47:49
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (35->20), mult. (83->34), div. (0->0), fcn. (91->6), ass. (0->18)
	t16 = pkin(8) + r_i_i_C(3);
	t6 = sin(qJ(2));
	t14 = t16 * t6;
	t9 = cos(qJ(2));
	t18 = t9 * pkin(2) + pkin(1) + t14;
	t7 = sin(qJ(1));
	t17 = t7 * t9;
	t10 = cos(qJ(1));
	t15 = t10 * t9;
	t5 = sin(qJ(3));
	t8 = cos(qJ(3));
	t12 = r_i_i_C(1) * t8 - r_i_i_C(2) * t5 + pkin(2);
	t11 = -t12 * t6 + t16 * t9;
	t4 = t8 * t15 + t7 * t5;
	t3 = -t5 * t15 + t7 * t8;
	t2 = t10 * t5 - t8 * t17;
	t1 = t10 * t8 + t5 * t17;
	t13 = [t10 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t18 * t7, t11 * t10, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0, 0; t7 * pkin(7) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t18 * t10, t11 * t7, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0, 0; 0, t12 * t9 + t14, (-r_i_i_C(1) * t5 - r_i_i_C(2) * t8) * t6, 0, 0, 0;];
	Ja_transl = t13;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:47:48
	% EndTime: 2019-10-10 11:47:49
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (61->22), mult. (142->34), div. (0->0), fcn. (163->6), ass. (0->21)
	t11 = cos(qJ(2));
	t19 = pkin(8) + r_i_i_C(2);
	t8 = sin(qJ(2));
	t15 = t19 * t8;
	t23 = pkin(2) * t11 + pkin(1) + t15;
	t10 = cos(qJ(3));
	t16 = r_i_i_C(3) + qJ(4);
	t20 = pkin(3) + r_i_i_C(1);
	t7 = sin(qJ(3));
	t22 = t20 * t10 + t16 * t7 + pkin(2);
	t9 = sin(qJ(1));
	t21 = t9 * t7;
	t18 = t9 * t10;
	t12 = cos(qJ(1));
	t17 = t11 * t12;
	t13 = t19 * t11 - t22 * t8;
	t4 = t10 * t17 + t21;
	t3 = t7 * t17 - t18;
	t2 = t11 * t18 - t12 * t7;
	t1 = t10 * t12 + t11 * t21;
	t5 = [pkin(7) * t12 - t16 * t1 - t20 * t2 - t23 * t9, t13 * t12, t16 * t4 - t20 * t3, t3, 0, 0; t9 * pkin(7) + t23 * t12 + t16 * t3 + t20 * t4, t13 * t9, -t20 * t1 + t16 * t2, t1, 0, 0; 0, t22 * t11 + t15, (t16 * t10 - t20 * t7) * t8, t8 * t7, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:47:48
	% EndTime: 2019-10-10 11:47:49
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (119->39), mult. (291->61), div. (0->0), fcn. (354->8), ass. (0->30)
	t18 = cos(qJ(2));
	t14 = sin(qJ(2));
	t28 = -r_i_i_C(3) - pkin(9) + pkin(8);
	t25 = t28 * t14;
	t34 = t18 * pkin(2) + pkin(1) + t25;
	t13 = sin(qJ(3));
	t17 = cos(qJ(3));
	t16 = cos(qJ(5));
	t12 = sin(qJ(5));
	t32 = pkin(3) + pkin(4);
	t24 = t12 * r_i_i_C(2) - t32;
	t21 = t16 * r_i_i_C(1) - t24;
	t26 = -t12 * r_i_i_C(1) - qJ(4);
	t22 = t16 * r_i_i_C(2) - t26;
	t33 = t22 * t13 + t21 * t17 + pkin(2);
	t19 = cos(qJ(1));
	t29 = t19 * t17;
	t15 = sin(qJ(1));
	t31 = t15 * t18;
	t6 = t13 * t31 + t29;
	t3 = t6 * t16;
	t30 = t19 * t13;
	t23 = (-(t12 * t17 - t13 * t16) * r_i_i_C(1) - (t12 * t13 + t16 * t17) * r_i_i_C(2)) * t14;
	t20 = -t33 * t14 + t28 * t18;
	t9 = t15 * t13 + t18 * t29;
	t8 = -t15 * t17 + t18 * t30;
	t7 = t17 * t31 - t30;
	t2 = t8 * t12 + t9 * t16;
	t1 = -t9 * t12 + t8 * t16;
	t4 = [t19 * pkin(7) - t3 * r_i_i_C(2) - t34 * t15 - t21 * t7 + t26 * t6, t20 * t19, -t21 * t8 + t22 * t9, t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t15 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t8 * qJ(4) + t34 * t19 + t32 * t9, t20 * t15, -t3 * r_i_i_C(1) + t22 * t7 + t24 * t6, t6, (-t7 * t12 + t3) * r_i_i_C(1) + (-t6 * t12 - t7 * t16) * r_i_i_C(2), 0; 0, t33 * t18 + t25, (qJ(4) * t17 - t32 * t13) * t14 - t23, t14 * t13, t23, 0;];
	Ja_transl = t4;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:47:48
	% EndTime: 2019-10-10 11:47:49
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (156->45), mult. (361->67), div. (0->0), fcn. (435->8), ass. (0->33)
	t15 = sin(qJ(3));
	t19 = cos(qJ(3));
	t18 = cos(qJ(5));
	t14 = sin(qJ(5));
	t37 = t18 * pkin(5) + pkin(3) + pkin(4);
	t29 = t14 * r_i_i_C(2) - t37;
	t24 = t18 * r_i_i_C(1) - t29;
	t39 = pkin(5) + r_i_i_C(1);
	t26 = -t14 * t39 - qJ(4);
	t41 = -t18 * r_i_i_C(2) + t26;
	t44 = t41 * t15 - t24 * t19 - pkin(2);
	t20 = cos(qJ(2));
	t16 = sin(qJ(2));
	t33 = -r_i_i_C(3) - qJ(6) - pkin(9) + pkin(8);
	t30 = t33 * t16;
	t43 = t20 * pkin(2) + pkin(1) + t30;
	t42 = (t14 * t19 - t15 * t18) * t16;
	t21 = cos(qJ(1));
	t34 = t21 * t19;
	t17 = sin(qJ(1));
	t36 = t17 * t20;
	t6 = t15 * t36 + t34;
	t3 = t6 * t18;
	t35 = t21 * t15;
	t31 = t14 * pkin(5) + qJ(4);
	t28 = -r_i_i_C(1) * t42 - (t14 * t15 + t18 * t19) * t16 * r_i_i_C(2);
	t8 = -t17 * t19 + t20 * t35;
	t9 = t17 * t15 + t20 * t34;
	t1 = -t9 * t14 + t8 * t18;
	t22 = t44 * t16 + t33 * t20;
	t7 = t19 * t36 - t35;
	t2 = t8 * t14 + t9 * t18;
	t4 = [t21 * pkin(7) - t3 * r_i_i_C(2) - t43 * t17 - t24 * t7 + t26 * t6, t22 * t21, -t24 * t8 - t41 * t9, t8, -t2 * r_i_i_C(2) + t39 * t1, -t21 * t16; t17 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t43 * t21 + t31 * t8 + t37 * t9, t22 * t17, -t3 * r_i_i_C(1) + t29 * t6 - t41 * t7, t6, (-t6 * t14 - t7 * t18) * r_i_i_C(2) + t39 * (-t7 * t14 + t3), -t17 * t16; 0, -t44 * t20 + t30, (-t15 * t37 + t19 * t31) * t16 - t28, t16 * t15, -pkin(5) * t42 + t28, t20;];
	Ja_transl = t4;
else
	Ja_transl=NaN(3,6);
end