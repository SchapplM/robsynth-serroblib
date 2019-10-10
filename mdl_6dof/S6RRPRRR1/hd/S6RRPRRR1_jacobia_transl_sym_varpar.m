% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:52
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:52:18
	% EndTime: 2019-10-10 10:52:18
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:52:18
	% EndTime: 2019-10-10 10:52:18
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
	% StartTime: 2019-10-10 10:52:18
	% EndTime: 2019-10-10 10:52:18
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
	% StartTime: 2019-10-10 10:52:18
	% EndTime: 2019-10-10 10:52:18
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (29->10), mult. (31->12), div. (0->0), fcn. (33->6), ass. (0->10)
	t5 = qJ(2) + pkin(11);
	t2 = sin(t5);
	t3 = cos(t5);
	t14 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2) + cos(qJ(2)) * pkin(2);
	t13 = r_i_i_C(3) + qJ(3) + pkin(7);
	t11 = pkin(1) + t14;
	t10 = -sin(qJ(2)) * pkin(2) - r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t11 * t8 + t13 * t9, t10 * t9, t8, 0, 0, 0; t11 * t9 + t13 * t8, t10 * t8, -t9, 0, 0, 0; 0, t14, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:52:18
	% EndTime: 2019-10-10 10:52:18
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (66->13), mult. (46->16), div. (0->0), fcn. (48->8), ass. (0->13)
	t10 = qJ(2) + pkin(11);
	t7 = qJ(4) + t10;
	t5 = sin(t7);
	t6 = cos(t7);
	t16 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t20 = t16 + pkin(3) * cos(t10) + cos(qJ(2)) * pkin(2);
	t18 = r_i_i_C(3) + pkin(8) + qJ(3) + pkin(7);
	t15 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
	t14 = pkin(1) + t20;
	t13 = t15 - pkin(3) * sin(t10) - sin(qJ(2)) * pkin(2);
	t12 = cos(qJ(1));
	t11 = sin(qJ(1));
	t1 = [-t14 * t11 + t18 * t12, t13 * t12, t11, t15 * t12, 0, 0; t18 * t11 + t14 * t12, t13 * t11, -t12, t15 * t11, 0, 0; 0, t20, 0, t16, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:52:18
	% EndTime: 2019-10-10 10:52:18
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (129->16), mult. (66->20), div. (0->0), fcn. (68->10), ass. (0->16)
	t13 = qJ(2) + pkin(11);
	t10 = qJ(4) + t13;
	t9 = qJ(5) + t10;
	t5 = sin(t9);
	t6 = cos(t9);
	t21 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t20 = t21 + pkin(4) * cos(t10);
	t26 = t20 + cos(qJ(2)) * pkin(2) + pkin(3) * cos(t13);
	t19 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
	t16 = t19 - pkin(4) * sin(t10);
	t22 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(3) + pkin(7);
	t18 = pkin(1) + t26;
	t17 = -pkin(3) * sin(t13) - sin(qJ(2)) * pkin(2) + t16;
	t15 = cos(qJ(1));
	t14 = sin(qJ(1));
	t1 = [-t18 * t14 + t22 * t15, t17 * t15, t14, t16 * t15, t19 * t15, 0; t22 * t14 + t18 * t15, t17 * t14, -t15, t16 * t14, t19 * t14, 0; 0, t26, 0, t20, t21, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:52:18
	% EndTime: 2019-10-10 10:52:18
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (288->37), mult. (167->45), div. (0->0), fcn. (177->12), ass. (0->34)
	t26 = qJ(2) + pkin(11);
	t23 = qJ(4) + t26;
	t22 = qJ(5) + t23;
	t18 = sin(t22);
	t19 = cos(t22);
	t27 = sin(qJ(6));
	t46 = r_i_i_C(2) * t27;
	t54 = pkin(10) + r_i_i_C(3);
	t55 = t18 * t46 + t19 * t54;
	t52 = t19 * pkin(5) + t54 * t18;
	t51 = pkin(3) * cos(t26) + cos(qJ(2)) * pkin(2);
	t29 = cos(qJ(6));
	t47 = r_i_i_C(1) * t29;
	t35 = (-pkin(5) - t47) * t18;
	t32 = t35 - pkin(4) * sin(t23);
	t17 = pkin(4) * cos(t23);
	t50 = pkin(1) + t17 + t51 + t52;
	t30 = cos(qJ(1));
	t43 = t27 * t30;
	t28 = sin(qJ(1));
	t42 = t28 * t27;
	t41 = t28 * t29;
	t40 = t29 * t30;
	t39 = t55 * t28;
	t37 = t55 * t30;
	t34 = -pkin(3) * sin(t26) - sin(qJ(2)) * pkin(2) + t32;
	t33 = (-t46 + t47) * t19 + t52;
	t31 = t17 + t33;
	t24 = -pkin(9) - pkin(8) - qJ(3) - pkin(7);
	t5 = t19 * t40 + t42;
	t4 = -t19 * t43 + t41;
	t3 = -t19 * t41 + t43;
	t2 = t19 * t42 + t40;
	t1 = [t3 * r_i_i_C(1) + t2 * r_i_i_C(2) - t24 * t30 - t50 * t28, t30 * t34 + t37, t28, t32 * t30 + t37, t30 * t35 + t37, r_i_i_C(1) * t4 - r_i_i_C(2) * t5; t5 * r_i_i_C(1) + t4 * r_i_i_C(2) - t28 * t24 + t50 * t30, t28 * t34 + t39, -t30, t28 * t32 + t39, t28 * t35 + t39, -r_i_i_C(1) * t2 + r_i_i_C(2) * t3; 0, t31 + t51, 0, t31, t33, (-r_i_i_C(1) * t27 - r_i_i_C(2) * t29) * t18;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end