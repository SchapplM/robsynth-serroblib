% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
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
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (9->6), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(11)) - r_i_i_C(2) * sin(pkin(11)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t5 * t3 + t6 * t4, t3, 0, 0, 0, 0; t6 * t3 + t5 * t4, -t4, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (26->9), mult. (26->11), div. (0->0), fcn. (28->5), ass. (0->10)
	t11 = r_i_i_C(3) + pkin(7) + qJ(2);
	t4 = pkin(11) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t10 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t9 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t8 = cos(pkin(11)) * pkin(2) + pkin(1) + t10;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [t11 * t7 - t8 * t6, t6, t9 * t7, 0, 0, 0; t11 * t6 + t8 * t7, -t7, t9 * t6, 0, 0, 0; 0, 0, t10, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (63->12), mult. (43->15), div. (0->0), fcn. (45->7), ass. (0->13)
	t9 = pkin(11) + qJ(3);
	t7 = qJ(4) + t9;
	t4 = sin(t7);
	t5 = cos(t7);
	t15 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t18 = t15 + pkin(3) * cos(t9);
	t16 = r_i_i_C(3) + pkin(8) + pkin(7) + qJ(2);
	t14 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5;
	t13 = cos(pkin(11)) * pkin(2) + pkin(1) + t18;
	t12 = -pkin(3) * sin(t9) + t14;
	t11 = cos(qJ(1));
	t10 = sin(qJ(1));
	t1 = [-t10 * t13 + t11 * t16, t10, t12 * t11, t14 * t11, 0, 0; t10 * t16 + t11 * t13, -t11, t12 * t10, t14 * t10, 0, 0; 0, 0, t18, t15, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (126->15), mult. (63->19), div. (0->0), fcn. (65->9), ass. (0->16)
	t12 = pkin(11) + qJ(3);
	t10 = qJ(4) + t12;
	t9 = qJ(5) + t10;
	t5 = sin(t9);
	t6 = cos(t9);
	t20 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t19 = t20 + pkin(4) * cos(t10);
	t24 = t19 + pkin(3) * cos(t12);
	t18 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
	t15 = t18 - pkin(4) * sin(t10);
	t21 = r_i_i_C(3) + pkin(9) + pkin(8) + pkin(7) + qJ(2);
	t17 = cos(pkin(11)) * pkin(2) + pkin(1) + t24;
	t16 = -pkin(3) * sin(t12) + t15;
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [-t17 * t13 + t21 * t14, t13, t16 * t14, t15 * t14, t18 * t14, 0; t21 * t13 + t17 * t14, -t14, t16 * t13, t15 * t13, t18 * t13, 0; 0, 0, t24, t19, t20, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (285->38), mult. (164->44), div. (0->0), fcn. (174->11), ass. (0->34)
	t25 = pkin(11) + qJ(3);
	t23 = qJ(4) + t25;
	t22 = qJ(5) + t23;
	t18 = sin(t22);
	t19 = cos(t22);
	t26 = sin(qJ(6));
	t45 = r_i_i_C(2) * t26;
	t52 = pkin(10) + r_i_i_C(3);
	t53 = t18 * t45 + t19 * t52;
	t50 = t19 * pkin(5) + t52 * t18;
	t28 = cos(qJ(6));
	t46 = r_i_i_C(1) * t28;
	t34 = (-pkin(5) - t46) * t18;
	t31 = t34 - pkin(4) * sin(t23);
	t17 = pkin(4) * cos(t23);
	t20 = pkin(3) * cos(t25);
	t49 = t17 + t20 + cos(pkin(11)) * pkin(2) + pkin(1) + t50;
	t29 = cos(qJ(1));
	t42 = t26 * t29;
	t27 = sin(qJ(1));
	t41 = t27 * t26;
	t40 = t27 * t28;
	t39 = t28 * t29;
	t38 = t53 * t27;
	t36 = t53 * t29;
	t33 = -pkin(3) * sin(t25) + t31;
	t32 = (-t45 + t46) * t19 + t50;
	t30 = t17 + t32;
	t24 = -pkin(9) - pkin(8) - pkin(7) - qJ(2);
	t5 = t19 * t39 + t41;
	t4 = -t19 * t42 + t40;
	t3 = -t19 * t40 + t42;
	t2 = t19 * t41 + t39;
	t1 = [t3 * r_i_i_C(1) + t2 * r_i_i_C(2) - t24 * t29 - t49 * t27, t27, t33 * t29 + t36, t31 * t29 + t36, t29 * t34 + t36, r_i_i_C(1) * t4 - r_i_i_C(2) * t5; t5 * r_i_i_C(1) + t4 * r_i_i_C(2) - t27 * t24 + t49 * t29, -t29, t33 * t27 + t38, t31 * t27 + t38, t27 * t34 + t38, -r_i_i_C(1) * t2 + r_i_i_C(2) * t3; 0, 0, t20 + t30, t30, t32, (-r_i_i_C(1) * t26 - r_i_i_C(2) * t28) * t18;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end