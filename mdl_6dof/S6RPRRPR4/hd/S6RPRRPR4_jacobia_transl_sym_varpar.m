% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:28
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
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
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (9->6), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(10)) - r_i_i_C(2) * sin(pkin(10)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t5 * t3 + t6 * t4, t3, 0, 0, 0, 0; t6 * t3 + t5 * t4, -t4, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (26->9), mult. (26->11), div. (0->0), fcn. (28->5), ass. (0->10)
	t11 = r_i_i_C(3) + pkin(7) + qJ(2);
	t4 = pkin(10) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t10 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t9 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t8 = cos(pkin(10)) * pkin(2) + pkin(1) + t10;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [t11 * t7 - t8 * t6, t6, t9 * t7, 0, 0, 0; t11 * t6 + t8 * t7, -t7, t9 * t6, 0, 0, 0; 0, 0, t10, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (63->12), mult. (43->15), div. (0->0), fcn. (45->7), ass. (0->13)
	t9 = pkin(10) + qJ(3);
	t7 = qJ(4) + t9;
	t4 = sin(t7);
	t5 = cos(t7);
	t15 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t18 = t15 + pkin(3) * cos(t9);
	t16 = r_i_i_C(3) + pkin(8) + pkin(7) + qJ(2);
	t14 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5;
	t13 = cos(pkin(10)) * pkin(2) + pkin(1) + t18;
	t12 = -pkin(3) * sin(t9) + t14;
	t11 = cos(qJ(1));
	t10 = sin(qJ(1));
	t1 = [-t10 * t13 + t11 * t16, t10, t12 * t11, t14 * t11, 0, 0; t10 * t16 + t11 * t13, -t11, t12 * t10, t14 * t10, 0, 0; 0, 0, t18, t15, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (144->23), mult. (105->24), div. (0->0), fcn. (114->9), ass. (0->21)
	t41 = r_i_i_C(3) + qJ(5);
	t18 = pkin(10) + qJ(3);
	t16 = qJ(4) + t18;
	t13 = sin(t16);
	t14 = cos(t16);
	t20 = cos(pkin(11));
	t28 = -r_i_i_C(1) * t20 - pkin(4);
	t19 = sin(pkin(11));
	t35 = r_i_i_C(2) * t19;
	t24 = t41 * t13 + (-t28 - t35) * t14;
	t40 = pkin(3) * cos(t18) + t24;
	t39 = t13 * t35 + t41 * t14;
	t37 = cos(pkin(10)) * pkin(2) + pkin(1) + t40;
	t21 = sin(qJ(1));
	t31 = t39 * t21;
	t22 = cos(qJ(1));
	t30 = t39 * t22;
	t27 = t28 * t13;
	t25 = t19 * r_i_i_C(1) + t20 * r_i_i_C(2) + pkin(7) + pkin(8) + qJ(2);
	t23 = -pkin(3) * sin(t18) + t27;
	t1 = [-t37 * t21 + t25 * t22, t21, t23 * t22 + t30, t22 * t27 + t30, t22 * t13, 0; t25 * t21 + t37 * t22, -t22, t23 * t21 + t31, t21 * t27 + t31, t21 * t13, 0; 0, 0, t40, t24, -t14, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (211->39), mult. (137->45), div. (0->0), fcn. (150->11), ass. (0->32)
	t23 = pkin(10) + qJ(3);
	t20 = qJ(4) + t23;
	t14 = sin(t20);
	t15 = cos(t20);
	t22 = pkin(11) + qJ(6);
	t17 = sin(t22);
	t42 = r_i_i_C(2) * t17;
	t46 = r_i_i_C(3) * t15 + t14 * t42;
	t16 = cos(pkin(11)) * pkin(5) + pkin(4);
	t25 = -pkin(9) - qJ(5);
	t45 = t15 * t16 + (r_i_i_C(3) - t25) * t14;
	t13 = pkin(3) * cos(t23);
	t44 = t13 + cos(pkin(10)) * pkin(2) + pkin(1) + t45;
	t19 = cos(t22);
	t43 = r_i_i_C(1) * t19;
	t26 = sin(qJ(1));
	t39 = t46 * t26;
	t27 = cos(qJ(1));
	t38 = t46 * t27;
	t37 = t17 * t27;
	t36 = t19 * t27;
	t35 = t26 * t17;
	t34 = t26 * t19;
	t32 = pkin(5) * sin(pkin(11)) + pkin(7) + pkin(8) + qJ(2);
	t30 = -t15 * t25 + (-t16 - t43) * t14;
	t29 = (-t42 + t43) * t15 + t45;
	t28 = -pkin(3) * sin(t23) + t30;
	t4 = t15 * t36 + t35;
	t3 = -t15 * t37 + t34;
	t2 = -t15 * t34 + t37;
	t1 = t15 * t35 + t36;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t44 * t26 + t32 * t27, t26, t28 * t27 + t38, t30 * t27 + t38, t27 * t14, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t32 * t26 + t44 * t27, -t27, t28 * t26 + t39, t30 * t26 + t39, t26 * t14, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, 0, t13 + t29, t29, -t15, (-r_i_i_C(1) * t17 - r_i_i_C(2) * t19) * t14;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end