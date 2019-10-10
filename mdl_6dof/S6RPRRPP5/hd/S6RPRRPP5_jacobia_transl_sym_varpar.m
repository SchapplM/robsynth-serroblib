% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:16
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPP5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:16:40
	% EndTime: 2019-10-10 01:16:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:16:40
	% EndTime: 2019-10-10 01:16:40
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
	% StartTime: 2019-10-10 01:16:40
	% EndTime: 2019-10-10 01:16:40
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (9->6), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(9)) - r_i_i_C(2) * sin(pkin(9)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t5 * t3 + t6 * t4, t3, 0, 0, 0, 0; t6 * t3 + t5 * t4, -t4, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:16:40
	% EndTime: 2019-10-10 01:16:40
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (26->9), mult. (26->11), div. (0->0), fcn. (28->5), ass. (0->10)
	t11 = r_i_i_C(3) + pkin(7) + qJ(2);
	t4 = pkin(9) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t10 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t9 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t8 = cos(pkin(9)) * pkin(2) + pkin(1) + t10;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [t11 * t7 - t8 * t6, t6, t9 * t7, 0, 0, 0; t11 * t6 + t8 * t7, -t7, t9 * t6, 0, 0, 0; 0, 0, t10, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:16:40
	% EndTime: 2019-10-10 01:16:40
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (71->24), mult. (85->33), div. (0->0), fcn. (95->7), ass. (0->22)
	t22 = pkin(8) + r_i_i_C(3);
	t8 = pkin(9) + qJ(3);
	t6 = sin(t8);
	t17 = t22 * t6;
	t7 = cos(t8);
	t23 = t17 + t7 * pkin(3) + cos(pkin(9)) * pkin(2) + pkin(1);
	t10 = sin(qJ(4));
	t13 = cos(qJ(1));
	t21 = t10 * t13;
	t11 = sin(qJ(1));
	t20 = t11 * t10;
	t12 = cos(qJ(4));
	t19 = t11 * t12;
	t18 = t12 * t13;
	t15 = r_i_i_C(1) * t12 - r_i_i_C(2) * t10 + pkin(3);
	t14 = -t15 * t6 + t22 * t7;
	t9 = -pkin(7) - qJ(2);
	t4 = t7 * t18 + t20;
	t3 = -t7 * t21 + t19;
	t2 = -t7 * t19 + t21;
	t1 = t7 * t20 + t18;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t23 * t11 - t13 * t9, t11, t14 * t13, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t11 * t9 + t23 * t13, -t13, t14 * t11, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, 0, t15 * t7 + t17, (-r_i_i_C(1) * t10 - r_i_i_C(2) * t12) * t6, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:16:40
	% EndTime: 2019-10-10 01:16:40
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (116->26), mult. (144->34), div. (0->0), fcn. (167->7), ass. (0->24)
	t24 = pkin(8) + r_i_i_C(2);
	t10 = pkin(9) + qJ(3);
	t8 = sin(t10);
	t18 = t24 * t8;
	t9 = cos(t10);
	t27 = t18 + pkin(3) * t9 + cos(pkin(9)) * pkin(2) + pkin(1);
	t12 = sin(qJ(4));
	t14 = cos(qJ(4));
	t19 = r_i_i_C(3) + qJ(5);
	t25 = pkin(4) + r_i_i_C(1);
	t26 = t19 * t12 + t25 * t14 + pkin(3);
	t13 = sin(qJ(1));
	t23 = t13 * t12;
	t22 = t13 * t14;
	t15 = cos(qJ(1));
	t21 = t14 * t15;
	t20 = t15 * t12;
	t16 = t24 * t9 - t26 * t8;
	t11 = -pkin(7) - qJ(2);
	t4 = t9 * t21 + t23;
	t3 = t9 * t20 - t22;
	t2 = t9 * t22 - t20;
	t1 = t9 * t23 + t21;
	t5 = [-t19 * t1 - t11 * t15 - t27 * t13 - t25 * t2, t13, t16 * t15, t19 * t4 - t25 * t3, t3, 0; -t13 * t11 + t27 * t15 + t19 * t3 + t25 * t4, -t15, t16 * t13, -t25 * t1 + t19 * t2, t1, 0; 0, 0, t26 * t9 + t18, (-t25 * t12 + t19 * t14) * t8, t8 * t12, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:16:40
	% EndTime: 2019-10-10 01:16:40
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (151->28), mult. (181->36), div. (0->0), fcn. (211->7), ass. (0->24)
	t19 = pkin(8) - r_i_i_C(3) - qJ(6);
	t10 = pkin(9) + qJ(3);
	t8 = sin(t10);
	t17 = t19 * t8;
	t9 = cos(t10);
	t27 = t17 + pkin(3) * t9 + cos(pkin(9)) * pkin(2) + pkin(1);
	t12 = sin(qJ(4));
	t14 = cos(qJ(4));
	t20 = pkin(4) + pkin(5) + r_i_i_C(1);
	t21 = r_i_i_C(2) + qJ(5);
	t26 = t21 * t12 + t20 * t14 + pkin(3);
	t13 = sin(qJ(1));
	t25 = t13 * t12;
	t24 = t13 * t14;
	t15 = cos(qJ(1));
	t23 = t14 * t15;
	t22 = t15 * t12;
	t16 = t19 * t9 - t26 * t8;
	t11 = -pkin(7) - qJ(2);
	t4 = t9 * t23 + t25;
	t3 = t9 * t22 - t24;
	t2 = t9 * t24 - t22;
	t1 = t9 * t25 + t23;
	t5 = [-t21 * t1 - t11 * t15 - t27 * t13 - t20 * t2, t13, t16 * t15, -t20 * t3 + t21 * t4, t3, -t15 * t8; -t13 * t11 + t27 * t15 + t20 * t4 + t21 * t3, -t15, t16 * t13, -t20 * t1 + t21 * t2, t1, -t13 * t8; 0, 0, t26 * t9 + t17, (-t20 * t12 + t21 * t14) * t8, t8 * t12, t9;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end