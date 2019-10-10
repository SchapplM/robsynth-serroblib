% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPP2
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
% Datum: 2019-10-10 01:11
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPP2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
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
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (8->5), mult. (6->6), div. (0->0), fcn. (6->4), ass. (0->4)
	t3 = qJ(1) + pkin(9);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [-t1 * r_i_i_C(1) - t2 * r_i_i_C(2) - sin(qJ(1)) * pkin(1), 0, 0, 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + cos(qJ(1)) * pkin(1), 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (27->9), mult. (26->12), div. (0->0), fcn. (26->6), ass. (0->10)
	t9 = pkin(7) + r_i_i_C(3);
	t4 = sin(qJ(3));
	t5 = cos(qJ(3));
	t8 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t7 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5;
	t6 = pkin(2) + t8;
	t3 = qJ(1) + pkin(9);
	t2 = cos(t3);
	t1 = sin(t3);
	t10 = [-sin(qJ(1)) * pkin(1) + t9 * t2 - t6 * t1, 0, t7 * t2, 0, 0, 0; cos(qJ(1)) * pkin(1) + t9 * t1 + t6 * t2, 0, t7 * t1, 0, 0, 0; 0, 1, t8, 0, 0, 0;];
	Ja_transl = t10;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (73->23), mult. (85->36), div. (0->0), fcn. (93->8), ass. (0->19)
	t11 = cos(qJ(3));
	t18 = pkin(8) + r_i_i_C(3);
	t9 = sin(qJ(3));
	t15 = t18 * t9;
	t19 = t11 * pkin(3) + pkin(2) + t15;
	t8 = sin(qJ(4));
	t17 = t11 * t8;
	t10 = cos(qJ(4));
	t16 = t10 * t11;
	t13 = r_i_i_C(1) * t10 - r_i_i_C(2) * t8 + pkin(3);
	t12 = t18 * t11 - t13 * t9;
	t7 = qJ(1) + pkin(9);
	t6 = cos(t7);
	t5 = sin(t7);
	t4 = t6 * t16 + t5 * t8;
	t3 = t5 * t10 - t6 * t17;
	t2 = -t5 * t16 + t6 * t8;
	t1 = t6 * t10 + t5 * t17;
	t14 = [-sin(qJ(1)) * pkin(1) + t6 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t19 * t5, 0, t12 * t6, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0; cos(qJ(1)) * pkin(1) + t5 * pkin(7) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t19 * t6, 0, t12 * t5, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, 1, t13 * t11 + t15, (-r_i_i_C(1) * t8 - r_i_i_C(2) * t10) * t9, 0, 0;];
	Ja_transl = t14;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (123->25), mult. (144->37), div. (0->0), fcn. (165->8), ass. (0->21)
	t13 = cos(qJ(3));
	t11 = sin(qJ(3));
	t20 = pkin(8) + r_i_i_C(2);
	t16 = t20 * t11;
	t23 = pkin(3) * t13 + pkin(2) + t16;
	t10 = sin(qJ(4));
	t12 = cos(qJ(4));
	t17 = r_i_i_C(3) + qJ(5);
	t21 = pkin(4) + r_i_i_C(1);
	t22 = t17 * t10 + t21 * t12 + pkin(3);
	t19 = t10 * t13;
	t18 = t12 * t13;
	t14 = -t22 * t11 + t20 * t13;
	t9 = qJ(1) + pkin(9);
	t8 = cos(t9);
	t7 = sin(t9);
	t4 = t7 * t10 + t8 * t18;
	t3 = -t7 * t12 + t8 * t19;
	t2 = -t8 * t10 + t7 * t18;
	t1 = t8 * t12 + t7 * t19;
	t5 = [-sin(qJ(1)) * pkin(1) + t8 * pkin(7) - t21 * t2 - t17 * t1 - t23 * t7, 0, t14 * t8, t17 * t4 - t21 * t3, t3, 0; cos(qJ(1)) * pkin(1) + t7 * pkin(7) + t21 * t4 + t17 * t3 + t23 * t8, 0, t14 * t7, -t21 * t1 + t17 * t2, t1, 0; 0, 1, t22 * t13 + t16, (-t21 * t10 + t17 * t12) * t11, t11 * t10, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (158->27), mult. (181->39), div. (0->0), fcn. (209->8), ass. (0->21)
	t13 = cos(qJ(3));
	t11 = sin(qJ(3));
	t17 = pkin(8) - r_i_i_C(3) - qJ(6);
	t15 = t17 * t11;
	t23 = t13 * pkin(3) + pkin(2) + t15;
	t10 = sin(qJ(4));
	t12 = cos(qJ(4));
	t18 = pkin(4) + pkin(5) + r_i_i_C(1);
	t19 = r_i_i_C(2) + qJ(5);
	t22 = t19 * t10 + t18 * t12 + pkin(3);
	t21 = t10 * t13;
	t20 = t12 * t13;
	t14 = -t22 * t11 + t17 * t13;
	t9 = qJ(1) + pkin(9);
	t8 = cos(t9);
	t7 = sin(t9);
	t4 = t7 * t10 + t8 * t20;
	t3 = -t7 * t12 + t8 * t21;
	t2 = -t8 * t10 + t7 * t20;
	t1 = t8 * t12 + t7 * t21;
	t5 = [-sin(qJ(1)) * pkin(1) + t8 * pkin(7) - t19 * t1 - t18 * t2 - t23 * t7, 0, t14 * t8, -t18 * t3 + t19 * t4, t3, -t8 * t11; cos(qJ(1)) * pkin(1) + t7 * pkin(7) + t19 * t3 + t18 * t4 + t23 * t8, 0, t14 * t7, -t18 * t1 + t19 * t2, t1, -t7 * t11; 0, 1, t22 * t13 + t15, (-t18 * t10 + t19 * t12) * t11, t11 * t10, t13;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end