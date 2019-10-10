% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:59
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRP8_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP8_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP8_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:59:40
	% EndTime: 2019-10-09 23:59:41
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:59:40
	% EndTime: 2019-10-09 23:59:41
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
	% StartTime: 2019-10-09 23:59:40
	% EndTime: 2019-10-09 23:59:40
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (7->5), mult. (8->4), div. (0->0), fcn. (10->2), ass. (0->5)
	t4 = pkin(1) - r_i_i_C(2);
	t3 = r_i_i_C(3) + qJ(2);
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t5 = [-t1 * t4 + t2 * t3, t1, 0, 0, 0, 0; t1 * t3 + t2 * t4, -t2, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:59:40
	% EndTime: 2019-10-09 23:59:41
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (11->6), mult. (16->6), div. (0->0), fcn. (20->4), ass. (0->5)
	t6 = pkin(1) + r_i_i_C(3) + qJ(3);
	t5 = r_i_i_C(1) * sin(pkin(9)) + r_i_i_C(2) * cos(pkin(9)) + qJ(2);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t6 * t3 + t5 * t4, t3, t4, 0, 0, 0; t5 * t3 + t6 * t4, -t4, t3, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:59:40
	% EndTime: 2019-10-09 23:59:41
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (28->10), mult. (30->11), div. (0->0), fcn. (34->5), ass. (0->10)
	t11 = pkin(1) + r_i_i_C(3) + pkin(7) + qJ(3);
	t3 = pkin(9) + qJ(4);
	t1 = sin(t3);
	t2 = cos(t3);
	t10 = r_i_i_C(1) * t2 - r_i_i_C(2) * t1;
	t9 = -t1 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t8 = pkin(3) * sin(pkin(9)) + qJ(2) - t9;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t4 = [-t11 * t6 + t8 * t7, t6, t7, t10 * t6, 0, 0; t11 * t7 + t8 * t6, -t7, t6, -t10 * t7, 0, 0; 0, 0, 0, t9, 0, 0;];
	Ja_transl = t4;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:59:40
	% EndTime: 2019-10-09 23:59:41
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (73->25), mult. (89->33), div. (0->0), fcn. (101->7), ass. (0->22)
	t20 = pkin(8) + r_i_i_C(3);
	t7 = pkin(9) + qJ(4);
	t6 = cos(t7);
	t23 = t20 * t6;
	t10 = sin(qJ(5));
	t12 = cos(qJ(5));
	t15 = r_i_i_C(1) * t12 - r_i_i_C(2) * t10 + pkin(4);
	t5 = sin(t7);
	t22 = t15 * t6 + t20 * t5;
	t21 = pkin(1) + pkin(7) + qJ(3);
	t13 = cos(qJ(1));
	t19 = t10 * t13;
	t11 = sin(qJ(1));
	t18 = t11 * t10;
	t17 = t11 * t12;
	t16 = t12 * t13;
	t14 = pkin(3) * sin(pkin(9)) + t5 * pkin(4) - t23 + qJ(2);
	t4 = t5 * t16 - t18;
	t3 = t5 * t19 + t17;
	t2 = t5 * t17 + t19;
	t1 = -t5 * t18 + t16;
	t8 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t21 * t11 + t14 * t13, t11, t13, t22 * t11, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t14 * t11 + t21 * t13, -t13, t11, -t22 * t13, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0; 0, 0, 0, -t15 * t5 + t23, (-r_i_i_C(1) * t10 - r_i_i_C(2) * t12) * t6, 0;];
	Ja_transl = t8;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:59:40
	% EndTime: 2019-10-09 23:59:41
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (118->28), mult. (148->34), div. (0->0), fcn. (173->7), ass. (0->24)
	t23 = pkin(8) + r_i_i_C(2);
	t11 = sin(qJ(5));
	t13 = cos(qJ(5));
	t17 = r_i_i_C(3) + qJ(6);
	t24 = pkin(5) + r_i_i_C(1);
	t25 = t17 * t11 + t24 * t13 + pkin(4);
	t8 = pkin(9) + qJ(4);
	t6 = sin(t8);
	t7 = cos(t8);
	t28 = t23 * t6 + t25 * t7;
	t26 = t23 * t7;
	t22 = pkin(1) + pkin(7) + qJ(3);
	t14 = cos(qJ(1));
	t21 = t11 * t14;
	t12 = sin(qJ(1));
	t20 = t12 * t11;
	t19 = t12 * t13;
	t18 = t14 * t13;
	t15 = pkin(3) * sin(pkin(9)) + pkin(4) * t6 - t26 + qJ(2);
	t4 = t6 * t18 - t20;
	t3 = t6 * t21 + t19;
	t2 = t6 * t19 + t21;
	t1 = t6 * t20 - t18;
	t5 = [-t22 * t12 + t15 * t14 + t17 * t3 + t24 * t4, t12, t14, t28 * t12, -t24 * t1 + t17 * t2, t1; t17 * t1 + t15 * t12 + t22 * t14 + t24 * t2, -t14, t12, -t28 * t14, -t17 * t4 + t24 * t3, -t3; 0, 0, 0, -t25 * t6 + t26, (-t24 * t11 + t17 * t13) * t7, t7 * t11;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end