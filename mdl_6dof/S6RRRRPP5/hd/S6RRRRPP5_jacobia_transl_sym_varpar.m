% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:27
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPP5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:27:32
	% EndTime: 2019-10-10 12:27:32
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:27:32
	% EndTime: 2019-10-10 12:27:32
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
	% StartTime: 2019-10-10 12:27:32
	% EndTime: 2019-10-10 12:27:32
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
	% StartTime: 2019-10-10 12:27:32
	% EndTime: 2019-10-10 12:27:32
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
	% StartTime: 2019-10-10 12:27:32
	% EndTime: 2019-10-10 12:27:32
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
	% StartTime: 2019-10-10 12:27:32
	% EndTime: 2019-10-10 12:27:32
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (168->31), mult. (203->44), div. (0->0), fcn. (232->8), ass. (0->31)
	t41 = pkin(4) + r_i_i_C(1);
	t34 = r_i_i_C(3) + qJ(5);
	t22 = cos(qJ(3));
	t15 = t22 * pkin(3) + pkin(2);
	t23 = cos(qJ(2));
	t20 = sin(qJ(2));
	t39 = r_i_i_C(2) + pkin(9) + pkin(8);
	t30 = t39 * t20;
	t43 = t23 * t15 + pkin(1) + t30;
	t18 = qJ(3) + qJ(4);
	t16 = sin(t18);
	t17 = cos(t18);
	t42 = t34 * t16 + t41 * t17 + t15;
	t19 = sin(qJ(3));
	t40 = pkin(3) * t19;
	t21 = sin(qJ(1));
	t37 = t21 * t23;
	t24 = cos(qJ(1));
	t36 = t24 * t16;
	t35 = t24 * t17;
	t33 = t34 * t17 * t20;
	t32 = t41 * t16;
	t31 = pkin(7) + t40;
	t7 = t16 * t37 + t35;
	t8 = t17 * t37 - t36;
	t28 = t34 * t8 - t41 * t7;
	t10 = t21 * t16 + t23 * t35;
	t9 = -t21 * t17 + t23 * t36;
	t27 = t34 * t10 - t41 * t9;
	t26 = -t42 * t20 + t39 * t23;
	t1 = [-t43 * t21 + t31 * t24 - t34 * t7 - t41 * t8, t26 * t24, (-t19 * t23 * t24 + t21 * t22) * pkin(3) + t27, t27, t9, 0; t41 * t10 + t31 * t21 + t43 * t24 + t34 * t9, t26 * t21, (-t19 * t37 - t22 * t24) * pkin(3) + t28, t28, t7, 0; 0, t42 * t23 + t30, (-t32 - t40) * t20 + t33, -t20 * t32 + t33, t20 * t16, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:27:32
	% EndTime: 2019-10-10 12:27:32
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (209->33), mult. (250->46), div. (0->0), fcn. (288->8), ass. (0->31)
	t38 = r_i_i_C(2) + qJ(5);
	t24 = cos(qJ(3));
	t17 = t24 * pkin(3) + pkin(2);
	t25 = cos(qJ(2));
	t22 = sin(qJ(2));
	t35 = r_i_i_C(3) + qJ(6) - pkin(9) - pkin(8);
	t31 = t35 * t22;
	t45 = -t25 * t17 - pkin(1) + t31;
	t36 = pkin(4) + pkin(5) + r_i_i_C(1);
	t20 = qJ(3) + qJ(4);
	t18 = sin(t20);
	t19 = cos(t20);
	t44 = t38 * t18 + t36 * t19 + t17;
	t21 = sin(qJ(3));
	t43 = pkin(3) * t21;
	t23 = sin(qJ(1));
	t41 = t23 * t25;
	t26 = cos(qJ(1));
	t40 = t26 * t18;
	t39 = t26 * t19;
	t37 = t38 * t19 * t22;
	t34 = pkin(7) + t43;
	t32 = t36 * t18;
	t10 = t19 * t41 - t40;
	t9 = t18 * t41 + t39;
	t30 = t38 * t10 - t36 * t9;
	t11 = -t23 * t19 + t25 * t40;
	t12 = t23 * t18 + t25 * t39;
	t29 = -t36 * t11 + t38 * t12;
	t28 = -t44 * t22 - t35 * t25;
	t1 = [-t36 * t10 + t45 * t23 + t34 * t26 - t38 * t9, t28 * t26, (-t21 * t25 * t26 + t23 * t24) * pkin(3) + t29, t29, t11, -t26 * t22; t38 * t11 + t36 * t12 + t34 * t23 - t45 * t26, t28 * t23, (-t21 * t41 - t24 * t26) * pkin(3) + t30, t30, t9, -t23 * t22; 0, t44 * t25 - t31, (-t32 - t43) * t22 + t37, -t22 * t32 + t37, t22 * t18, t25;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end