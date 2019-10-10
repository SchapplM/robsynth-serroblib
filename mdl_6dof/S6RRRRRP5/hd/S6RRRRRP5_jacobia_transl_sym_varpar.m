% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:03
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:03:55
	% EndTime: 2019-10-10 13:03:55
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:03:55
	% EndTime: 2019-10-10 13:03:55
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
	% StartTime: 2019-10-10 13:03:55
	% EndTime: 2019-10-10 13:03:55
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
	% StartTime: 2019-10-10 13:03:55
	% EndTime: 2019-10-10 13:03:55
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
	% StartTime: 2019-10-10 13:03:55
	% EndTime: 2019-10-10 13:03:55
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
	% StartTime: 2019-10-10 13:03:55
	% EndTime: 2019-10-10 13:03:55
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (193->36), mult. (168->50), div. (0->0), fcn. (184->10), ass. (0->30)
	t23 = cos(qJ(2));
	t21 = sin(qJ(2));
	t32 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8);
	t28 = t32 * t21;
	t20 = qJ(3) + qJ(4);
	t16 = cos(t20);
	t11 = pkin(4) * t16 + cos(qJ(3)) * pkin(3);
	t9 = pkin(2) + t11;
	t37 = t23 * t9 + pkin(1) + t28;
	t17 = qJ(5) + t20;
	t13 = sin(t17);
	t14 = cos(t17);
	t24 = cos(qJ(1));
	t22 = sin(qJ(1));
	t31 = t22 * t23;
	t5 = t13 * t31 + t14 * t24;
	t6 = t13 * t24 - t14 * t31;
	t36 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t30 = t23 * t24;
	t7 = -t13 * t30 + t22 * t14;
	t8 = t22 * t13 + t14 * t30;
	t35 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t15 = sin(t20);
	t34 = pkin(4) * t15;
	t10 = t34 + sin(qJ(3)) * pkin(3);
	t33 = pkin(7) + t10;
	t27 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
	t26 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t9;
	t25 = -t26 * t21 + t32 * t23;
	t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t37 * t22 + t33 * t24, t25 * t24, -t10 * t30 + t22 * t11 + t35, (-t15 * t30 + t16 * t22) * pkin(4) + t35, t35, 0; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t33 * t22 + t37 * t24, t25 * t22, -t10 * t31 - t11 * t24 + t36, (-t15 * t31 - t16 * t24) * pkin(4) + t36, t36, 0; 0, t26 * t23 + t28, (-t10 + t27) * t21, (t27 - t34) * t21, t27 * t21, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:03:55
	% EndTime: 2019-10-10 13:03:55
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (265->43), mult. (199->55), div. (0->0), fcn. (218->10), ass. (0->32)
	t25 = cos(qJ(2));
	t23 = sin(qJ(2));
	t36 = r_i_i_C(3) + qJ(6) + pkin(10) + pkin(9) + pkin(8);
	t30 = t36 * t23;
	t22 = qJ(3) + qJ(4);
	t20 = qJ(5) + t22;
	t17 = cos(t20);
	t13 = pkin(5) * t17 + pkin(4) * cos(t22);
	t11 = cos(qJ(3)) * pkin(3) + t13;
	t9 = pkin(2) + t11;
	t41 = t25 * t9 + pkin(1) + t30;
	t16 = sin(t20);
	t26 = cos(qJ(1));
	t32 = t26 * t17;
	t24 = sin(qJ(1));
	t35 = t24 * t25;
	t5 = t16 * t35 + t32;
	t33 = t26 * t16;
	t6 = -t17 * t35 + t33;
	t40 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t7 = t24 * t17 - t25 * t33;
	t8 = t24 * t16 + t25 * t32;
	t39 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t38 = r_i_i_C(2) * t17;
	t12 = -pkin(4) * sin(t22) - pkin(5) * t16;
	t10 = sin(qJ(3)) * pkin(3) - t12;
	t37 = pkin(7) + t10;
	t34 = t25 * t26;
	t29 = -r_i_i_C(1) * t16 - t38;
	t28 = r_i_i_C(1) * t17 - r_i_i_C(2) * t16 + t9;
	t27 = -t28 * t23 + t36 * t25;
	t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t41 * t24 + t37 * t26, t27 * t26, -t10 * t34 + t24 * t11 + t39, t12 * t34 + t24 * t13 + t39, t7 * pkin(5) + t39, t26 * t23; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t37 * t24 + t41 * t26, t27 * t24, -t10 * t35 - t26 * t11 + t40, t12 * t35 - t26 * t13 + t40, -t5 * pkin(5) + t40, t24 * t23; 0, t28 * t25 + t30, (-t10 + t29) * t23, (t12 + t29) * t23, (-t38 + (-pkin(5) - r_i_i_C(1)) * t16) * t23, -t25;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end