% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:04
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:25
	% EndTime: 2019-10-10 10:04:25
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:25
	% EndTime: 2019-10-10 10:04:25
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:25
	% EndTime: 2019-10-10 10:04:25
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
	% StartTime: 2019-10-10 10:04:25
	% EndTime: 2019-10-10 10:04:25
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (29->10), mult. (31->12), div. (0->0), fcn. (33->6), ass. (0->10)
	t5 = qJ(2) + pkin(10);
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
	% StartTime: 2019-10-10 10:04:25
	% EndTime: 2019-10-10 10:04:25
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (66->13), mult. (46->16), div. (0->0), fcn. (48->8), ass. (0->13)
	t10 = qJ(2) + pkin(10);
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
	% StartTime: 2019-10-10 10:04:25
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (147->24), mult. (108->25), div. (0->0), fcn. (117->10), ass. (0->21)
	t43 = r_i_i_C(3) + qJ(5);
	t19 = qJ(2) + pkin(10);
	t16 = qJ(4) + t19;
	t14 = sin(t16);
	t15 = cos(t16);
	t21 = cos(pkin(11));
	t29 = -r_i_i_C(1) * t21 - pkin(4);
	t20 = sin(pkin(11));
	t37 = r_i_i_C(2) * t20;
	t24 = t43 * t14 + (-t29 - t37) * t15;
	t42 = t24 + pkin(3) * cos(t19) + cos(qJ(2)) * pkin(2);
	t41 = t14 * t37 + t43 * t15;
	t39 = pkin(1) + t42;
	t22 = sin(qJ(1));
	t32 = t41 * t22;
	t23 = cos(qJ(1));
	t31 = t41 * t23;
	t28 = t29 * t14;
	t26 = t20 * r_i_i_C(1) + t21 * r_i_i_C(2) + pkin(7) + pkin(8) + qJ(3);
	t25 = -pkin(3) * sin(t19) - sin(qJ(2)) * pkin(2) + t28;
	t1 = [-t39 * t22 + t26 * t23, t25 * t23 + t31, t22, t23 * t28 + t31, t23 * t14, 0; t26 * t22 + t39 * t23, t25 * t22 + t32, -t23, t22 * t28 + t32, t22 * t14, 0; 0, t42, 0, t24, -t15, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:25
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (214->39), mult. (140->46), div. (0->0), fcn. (153->12), ass. (0->32)
	t24 = qJ(2) + pkin(10);
	t20 = qJ(4) + t24;
	t15 = sin(t20);
	t16 = cos(t20);
	t23 = pkin(11) + qJ(6);
	t18 = sin(t23);
	t44 = r_i_i_C(2) * t18;
	t48 = r_i_i_C(3) * t16 + t15 * t44;
	t17 = cos(pkin(11)) * pkin(5) + pkin(4);
	t26 = -pkin(9) - qJ(5);
	t47 = t16 * t17 + (r_i_i_C(3) - t26) * t15;
	t35 = pkin(3) * cos(t24) + cos(qJ(2)) * pkin(2);
	t46 = pkin(1) + t35 + t47;
	t19 = cos(t23);
	t45 = r_i_i_C(1) * t19;
	t27 = sin(qJ(1));
	t41 = t48 * t27;
	t28 = cos(qJ(1));
	t40 = t48 * t28;
	t39 = t18 * t28;
	t38 = t19 * t28;
	t37 = t27 * t18;
	t36 = t27 * t19;
	t33 = pkin(5) * sin(pkin(11)) + pkin(7) + pkin(8) + qJ(3);
	t31 = -t16 * t26 + (-t17 - t45) * t15;
	t30 = (-t44 + t45) * t16 + t47;
	t29 = -pkin(3) * sin(t24) - sin(qJ(2)) * pkin(2) + t31;
	t4 = t16 * t38 + t37;
	t3 = -t16 * t39 + t36;
	t2 = -t16 * t36 + t39;
	t1 = t16 * t37 + t38;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t46 * t27 + t33 * t28, t29 * t28 + t40, t27, t31 * t28 + t40, t28 * t15, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t33 * t27 + t46 * t28, t29 * t27 + t41, -t28, t31 * t27 + t41, t27 * t15, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, t30 + t35, 0, t30, -t16, (-r_i_i_C(1) * t18 - r_i_i_C(2) * t19) * t15;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end