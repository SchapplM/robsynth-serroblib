% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:25
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPP4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:45
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:45
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
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:45
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:45
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
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:45
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
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:45
	% DurationCPUTime: 0.24s
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
	t15 = pkin(10) + t20;
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
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:45
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (281->40), mult. (229->52), div. (0->0), fcn. (261->10), ass. (0->35)
	t47 = pkin(5) + r_i_i_C(1);
	t38 = r_i_i_C(3) + qJ(6);
	t26 = qJ(3) + qJ(4);
	t23 = cos(t26);
	t13 = pkin(4) * t23 + cos(qJ(3)) * pkin(3);
	t11 = pkin(2) + t13;
	t29 = cos(qJ(2));
	t27 = sin(qJ(2));
	t44 = r_i_i_C(2) + qJ(5) + pkin(9) + pkin(8);
	t35 = t44 * t27;
	t49 = t29 * t11 + pkin(1) + t35;
	t21 = pkin(10) + t26;
	t18 = sin(t21);
	t19 = cos(t21);
	t48 = t38 * t18 + t47 * t19 + t11;
	t22 = sin(t26);
	t46 = pkin(4) * t22;
	t12 = t46 + sin(qJ(3)) * pkin(3);
	t45 = pkin(7) + t12;
	t28 = sin(qJ(1));
	t42 = t28 * t29;
	t30 = cos(qJ(1));
	t41 = t29 * t30;
	t40 = t30 * t18;
	t39 = t30 * t19;
	t37 = t38 * t19 * t27;
	t36 = t47 * t18;
	t7 = t18 * t42 + t39;
	t8 = t19 * t42 - t40;
	t33 = t38 * t8 - t47 * t7;
	t10 = t28 * t18 + t29 * t39;
	t9 = -t28 * t19 + t29 * t40;
	t32 = t38 * t10 - t47 * t9;
	t31 = -t48 * t27 + t44 * t29;
	t1 = [-t49 * t28 + t45 * t30 - t38 * t7 - t47 * t8, t31 * t30, -t12 * t41 + t28 * t13 + t32, (-t22 * t41 + t23 * t28) * pkin(4) + t32, t30 * t27, t9; t47 * t10 + t45 * t28 + t49 * t30 + t38 * t9, t31 * t28, -t12 * t42 - t30 * t13 + t33, (-t22 * t42 - t23 * t30) * pkin(4) + t33, t28 * t27, t7; 0, t48 * t29 + t35, (-t12 - t36) * t27 + t37, (-t36 - t46) * t27 + t37, -t29, t27 * t18;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end