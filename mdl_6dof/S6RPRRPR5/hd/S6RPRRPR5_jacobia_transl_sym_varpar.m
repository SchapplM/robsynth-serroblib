% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
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
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
	% DurationCPUTime: 0.12s
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
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
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
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (118->20), mult. (73->21), div. (0->0), fcn. (78->7), ass. (0->18)
	t34 = r_i_i_C(3) + qJ(5);
	t17 = pkin(10) + qJ(3);
	t15 = qJ(4) + t17;
	t12 = sin(t15);
	t13 = cos(t15);
	t20 = (pkin(4) - r_i_i_C(2)) * t13 + t34 * t12;
	t33 = pkin(3) * cos(t17) + t20;
	t32 = t34 * t13;
	t31 = cos(pkin(10)) * pkin(2) + pkin(1) + t33;
	t28 = r_i_i_C(1) + pkin(8) + pkin(7) + qJ(2);
	t18 = sin(qJ(1));
	t27 = t18 * t12;
	t19 = cos(qJ(1));
	t26 = t19 * t12;
	t23 = r_i_i_C(2) * t27 + t32 * t18;
	t22 = r_i_i_C(2) * t26 + t32 * t19;
	t21 = -pkin(3) * sin(t17) - pkin(4) * t12;
	t1 = [-t31 * t18 + t28 * t19, t18, t21 * t19 + t22, -pkin(4) * t26 + t22, t26, 0; t28 * t18 + t31 * t19, -t19, t21 * t18 + t23, -pkin(4) * t27 + t23, t27, 0; 0, 0, t33, t20, -t13, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (189->32), mult. (141->41), div. (0->0), fcn. (154->9), ass. (0->28)
	t25 = sin(qJ(6));
	t27 = cos(qJ(6));
	t49 = r_i_i_C(1) * t25 + r_i_i_C(2) * t27;
	t24 = pkin(10) + qJ(3);
	t22 = qJ(4) + t24;
	t19 = sin(t22);
	t20 = cos(t22);
	t37 = pkin(4) + pkin(9) + r_i_i_C(3);
	t48 = t19 * qJ(5) + t37 * t20;
	t46 = (qJ(5) + t49) * t20;
	t18 = pkin(3) * cos(t24);
	t45 = t18 + cos(pkin(10)) * pkin(2) + pkin(1) + t48;
	t42 = pkin(5) + pkin(8) + pkin(7) + qJ(2);
	t26 = sin(qJ(1));
	t41 = t26 * t25;
	t40 = t26 * t27;
	t28 = cos(qJ(1));
	t39 = t28 * t19;
	t36 = t46 * t26;
	t35 = t46 * t28;
	t31 = t37 * t19;
	t30 = t49 * t19 + t48;
	t29 = -pkin(3) * sin(t24) - t31;
	t4 = -t19 * t41 + t27 * t28;
	t3 = t19 * t40 + t25 * t28;
	t2 = t25 * t39 + t40;
	t1 = t27 * t39 - t41;
	t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t45 * t26 + t42 * t28, t26, t29 * t28 + t35, -t28 * t31 + t35, t39, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t42 * t26 + t45 * t28, -t28, t29 * t26 + t36, -t26 * t31 + t36, t26 * t19, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 0, t18 + t30, t30, -t20, (-r_i_i_C(1) * t27 + r_i_i_C(2) * t25) * t20;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end