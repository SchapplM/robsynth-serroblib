% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:27
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR7_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR7_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR7_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:45
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:45
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
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:45
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
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:45
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
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:45
	% DurationCPUTime: 0.17s
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
	t5 = [pkin(7) * t12 - t16 * t1 - t20 * t2 - t23 * t9, t13 * t12, t16 * t4 - t20 * t3, t3, 0, 0; t9 * pkin(7) + t12 * t23 + t16 * t3 + t20 * t4, t13 * t9, -t20 * t1 + t16 * t2, t1, 0, 0; 0, t22 * t11 + t15, (t16 * t10 - t20 * t7) * t8, t8 * t7, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:45
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (104->31), mult. (247->46), div. (0->0), fcn. (299->8), ass. (0->27)
	t15 = cos(qJ(2));
	t12 = sin(qJ(2));
	t24 = -r_i_i_C(3) - qJ(5) + pkin(8);
	t20 = t24 * t12;
	t29 = t15 * pkin(2) + pkin(1) + t20;
	t11 = sin(qJ(3));
	t14 = cos(qJ(3));
	t10 = cos(pkin(10));
	t9 = sin(pkin(10));
	t21 = t9 * r_i_i_C(2) - pkin(3) - pkin(4);
	t18 = t10 * r_i_i_C(1) - t21;
	t22 = t9 * r_i_i_C(1) + qJ(4);
	t19 = t10 * r_i_i_C(2) + t22;
	t28 = t19 * t11 + t18 * t14 + pkin(2);
	t13 = sin(qJ(1));
	t27 = t13 * t15;
	t16 = cos(qJ(1));
	t26 = t16 * t11;
	t25 = t16 * t14;
	t17 = -t28 * t12 + t24 * t15;
	t6 = t13 * t11 + t15 * t25;
	t5 = -t13 * t14 + t15 * t26;
	t4 = t14 * t27 - t26;
	t3 = t11 * t27 + t25;
	t2 = t6 * t10;
	t1 = t3 * t10;
	t7 = [t16 * pkin(7) - t1 * r_i_i_C(2) - t29 * t13 - t18 * t4 - t22 * t3, t17 * t16, t2 * r_i_i_C(2) - t18 * t5 + t22 * t6, t5, -t16 * t12, 0; t13 * pkin(7) + t2 * r_i_i_C(1) + t29 * t16 + t19 * t5 - t21 * t6, t17 * t13, -t1 * r_i_i_C(1) + t19 * t4 + t21 * t3, t3, -t13 * t12, 0; 0, t28 * t15 + t20, (-t18 * t11 + t19 * t14) * t12, t12 * t11, t15, 0;];
	Ja_transl = t7;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:46
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (190->44), mult. (335->65), div. (0->0), fcn. (405->10), ass. (0->32)
	t22 = cos(qJ(2));
	t19 = sin(qJ(2));
	t33 = -r_i_i_C(3) - pkin(9) - qJ(5) + pkin(8);
	t30 = t33 * t19;
	t39 = t22 * pkin(2) + pkin(1) + t30;
	t18 = sin(qJ(3));
	t21 = cos(qJ(3));
	t15 = pkin(10) + qJ(6);
	t14 = cos(t15);
	t13 = sin(t15);
	t31 = sin(pkin(10)) * pkin(5) + qJ(4);
	t27 = -t13 * r_i_i_C(1) - t31;
	t25 = t14 * r_i_i_C(2) - t27;
	t37 = pkin(3) + cos(pkin(10)) * pkin(5) + pkin(4);
	t29 = t13 * r_i_i_C(2) - t37;
	t26 = t14 * r_i_i_C(1) - t29;
	t38 = t18 * t25 + t21 * t26 + pkin(2);
	t23 = cos(qJ(1));
	t34 = t23 * t21;
	t20 = sin(qJ(1));
	t36 = t20 * t22;
	t6 = t18 * t36 + t34;
	t5 = t6 * t14;
	t35 = t23 * t18;
	t28 = (-(t13 * t21 - t14 * t18) * r_i_i_C(1) - (t13 * t18 + t14 * t21) * r_i_i_C(2)) * t19;
	t24 = -t38 * t19 + t33 * t22;
	t9 = t20 * t18 + t22 * t34;
	t8 = -t20 * t21 + t22 * t35;
	t7 = t21 * t36 - t35;
	t2 = t8 * t13 + t9 * t14;
	t1 = -t9 * t13 + t8 * t14;
	t3 = [t23 * pkin(7) - t5 * r_i_i_C(2) - t39 * t20 - t26 * t7 + t27 * t6, t24 * t23, t25 * t9 - t26 * t8, t8, -t23 * t19, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t20 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t39 * t23 + t31 * t8 + t37 * t9, t24 * t20, -t5 * r_i_i_C(1) + t25 * t7 + t29 * t6, t6, -t20 * t19, (-t7 * t13 + t5) * r_i_i_C(1) + (-t6 * t13 - t7 * t14) * r_i_i_C(2); 0, t38 * t22 + t30, (-t37 * t18 + t31 * t21) * t19 - t28, t19 * t18, t22, t28;];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end