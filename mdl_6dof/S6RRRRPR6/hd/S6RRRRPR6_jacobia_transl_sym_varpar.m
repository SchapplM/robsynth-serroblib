% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:40
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:40:25
	% EndTime: 2019-10-10 12:40:25
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:40:25
	% EndTime: 2019-10-10 12:40:25
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
	% StartTime: 2019-10-10 12:40:25
	% EndTime: 2019-10-10 12:40:25
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (11->6), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->9)
	t8 = pkin(7) + r_i_i_C(3);
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t2 * t5 + t4 * t8, t6 * t4, 0, 0, 0, 0; t2 * t8 + t4 * t5, t6 * t2, 0, 0, 0, 0; 0, t7, 0, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:40:25
	% EndTime: 2019-10-10 12:40:25
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
	t13 = [t10 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t18 * t7, t11 * t10, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0, 0; t7 * pkin(7) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t18 * t10, t11 * t7, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0, 0; 0, t12 * t9 + t14, (-r_i_i_C(1) * t5 - r_i_i_C(2) * t8) * t6, 0, 0, 0;];
	Ja_transl = t13;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:40:25
	% EndTime: 2019-10-10 12:40:25
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (93->28), mult. (124->43), div. (0->0), fcn. (136->8), ass. (0->27)
	t17 = cos(qJ(2));
	t14 = sin(qJ(2));
	t28 = r_i_i_C(3) + pkin(9) + pkin(8);
	t23 = t28 * t14;
	t16 = cos(qJ(3));
	t9 = pkin(3) * t16 + pkin(2);
	t32 = t17 * t9 + pkin(1) + t23;
	t12 = qJ(3) + qJ(4);
	t10 = sin(t12);
	t11 = cos(t12);
	t18 = cos(qJ(1));
	t15 = sin(qJ(1));
	t27 = t15 * t17;
	t5 = t10 * t27 + t11 * t18;
	t6 = t10 * t18 - t11 * t27;
	t31 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t26 = t17 * t18;
	t7 = -t10 * t26 + t11 * t15;
	t8 = t10 * t15 + t11 * t26;
	t30 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t13 = sin(qJ(3));
	t29 = pkin(3) * t13;
	t24 = pkin(7) + t29;
	t22 = -r_i_i_C(1) * t10 - r_i_i_C(2) * t11;
	t21 = r_i_i_C(1) * t11 - r_i_i_C(2) * t10 + t9;
	t20 = -t21 * t14 + t28 * t17;
	t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t32 * t15 + t24 * t18, t20 * t18, (-t13 * t26 + t15 * t16) * pkin(3) + t30, t30, 0, 0; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t24 * t15 + t32 * t18, t20 * t15, (-t13 * t27 - t16 * t18) * pkin(3) + t31, t31, 0, 0; 0, t21 * t17 + t23, (t22 - t29) * t14, t22 * t14, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:40:25
	% EndTime: 2019-10-10 12:40:25
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (167->37), mult. (150->51), div. (0->0), fcn. (165->10), ass. (0->32)
	t23 = cos(qJ(2));
	t21 = sin(qJ(2));
	t34 = r_i_i_C(3) + qJ(5) + pkin(9) + pkin(8);
	t28 = t34 * t21;
	t20 = qJ(3) + qJ(4);
	t17 = cos(t20);
	t11 = pkin(4) * t17 + cos(qJ(3)) * pkin(3);
	t9 = pkin(2) + t11;
	t39 = t23 * t9 + pkin(1) + t28;
	t15 = pkin(11) + t20;
	t12 = sin(t15);
	t13 = cos(t15);
	t24 = cos(qJ(1));
	t30 = t24 * t13;
	t22 = sin(qJ(1));
	t33 = t22 * t23;
	t5 = t12 * t33 + t30;
	t31 = t24 * t12;
	t6 = -t13 * t33 + t31;
	t38 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t7 = t22 * t13 - t23 * t31;
	t8 = t22 * t12 + t23 * t30;
	t37 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t16 = sin(t20);
	t36 = pkin(4) * t16;
	t10 = t36 + sin(qJ(3)) * pkin(3);
	t35 = pkin(7) + t10;
	t32 = t23 * t24;
	t27 = -r_i_i_C(1) * t12 - r_i_i_C(2) * t13;
	t26 = r_i_i_C(1) * t13 - r_i_i_C(2) * t12 + t9;
	t25 = -t26 * t21 + t34 * t23;
	t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t39 * t22 + t35 * t24, t25 * t24, -t10 * t32 + t22 * t11 + t37, (-t16 * t32 + t17 * t22) * pkin(4) + t37, t24 * t21, 0; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t35 * t22 + t39 * t24, t25 * t22, -t10 * t33 - t24 * t11 + t38, (-t16 * t33 - t17 * t24) * pkin(4) + t38, t22 * t21, 0; 0, t26 * t23 + t28, (-t10 + t27) * t21, (t27 - t36) * t21, -t23, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:40:25
	% EndTime: 2019-10-10 12:40:25
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (294->40), mult. (187->52), div. (0->0), fcn. (206->12), ass. (0->32)
	t27 = cos(qJ(2));
	t25 = sin(qJ(2));
	t38 = r_i_i_C(3) + pkin(10) + qJ(5) + pkin(9) + pkin(8);
	t32 = t38 * t25;
	t24 = qJ(3) + qJ(4);
	t20 = pkin(11) + t24;
	t13 = pkin(5) * cos(t20) + pkin(4) * cos(t24);
	t11 = cos(qJ(3)) * pkin(3) + t13;
	t9 = pkin(2) + t11;
	t42 = t27 * t9 + pkin(1) + t32;
	t19 = qJ(6) + t20;
	t15 = sin(t19);
	t16 = cos(t19);
	t28 = cos(qJ(1));
	t34 = t28 * t16;
	t26 = sin(qJ(1));
	t37 = t26 * t27;
	t5 = t15 * t37 + t34;
	t35 = t28 * t15;
	t6 = -t16 * t37 + t35;
	t41 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t7 = t26 * t16 - t27 * t35;
	t8 = t26 * t15 + t27 * t34;
	t40 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t12 = -pkin(4) * sin(t24) - pkin(5) * sin(t20);
	t10 = sin(qJ(3)) * pkin(3) - t12;
	t39 = pkin(7) + t10;
	t36 = t27 * t28;
	t31 = -r_i_i_C(1) * t15 - r_i_i_C(2) * t16;
	t30 = r_i_i_C(1) * t16 - r_i_i_C(2) * t15 + t9;
	t29 = -t30 * t25 + t38 * t27;
	t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t42 * t26 + t39 * t28, t29 * t28, -t10 * t36 + t26 * t11 + t40, t12 * t36 + t26 * t13 + t40, t28 * t25, t40; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t39 * t26 + t42 * t28, t29 * t26, -t10 * t37 - t28 * t11 + t41, t12 * t37 - t28 * t13 + t41, t26 * t25, t41; 0, t30 * t27 + t32, (-t10 + t31) * t25, (t12 + t31) * t25, -t27, t31 * t25;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end