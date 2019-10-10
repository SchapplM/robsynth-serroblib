% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPP8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:21
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPP8_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP8_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP8_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:48
	% EndTime: 2019-10-10 01:21:48
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:48
	% EndTime: 2019-10-10 01:21:48
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
	% StartTime: 2019-10-10 01:21:48
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 01:21:48
	% EndTime: 2019-10-10 01:21:48
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (14->8), mult. (26->10), div. (0->0), fcn. (28->4), ass. (0->9)
	t8 = pkin(1) + pkin(7) + r_i_i_C(3);
	t1 = sin(qJ(3));
	t3 = cos(qJ(3));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -t1 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t5 = qJ(2) - t6;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t8 * t2 + t5 * t4, t2, t7 * t2, 0, 0, 0; t5 * t2 + t8 * t4, -t4, -t7 * t4, 0, 0, 0; 0, 0, t6, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:48
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (38->23), mult. (85->32), div. (0->0), fcn. (95->6), ass. (0->21)
	t15 = pkin(8) + r_i_i_C(3);
	t9 = cos(qJ(3));
	t20 = t15 * t9;
	t5 = sin(qJ(4));
	t8 = cos(qJ(4));
	t12 = r_i_i_C(1) * t8 - r_i_i_C(2) * t5 + pkin(3);
	t6 = sin(qJ(3));
	t19 = t12 * t9 + t15 * t6;
	t7 = sin(qJ(1));
	t18 = t7 * t5;
	t17 = t7 * t8;
	t16 = pkin(1) + pkin(7);
	t10 = cos(qJ(1));
	t14 = t10 * t5;
	t13 = t10 * t8;
	t11 = t6 * pkin(3) + qJ(2) - t20;
	t4 = t6 * t13 - t18;
	t3 = t6 * t14 + t17;
	t2 = t6 * t17 + t14;
	t1 = -t6 * t18 + t13;
	t21 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t11 * t10 - t16 * t7, t7, t19 * t7, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t16 * t10 + t11 * t7, -t10, -t19 * t10, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0, 0; 0, 0, -t12 * t6 + t20, (-r_i_i_C(1) * t5 - r_i_i_C(2) * t8) * t9, 0, 0;];
	Ja_transl = t21;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:48
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (64->26), mult. (144->33), div. (0->0), fcn. (167->6), ass. (0->23)
	t11 = cos(qJ(3));
	t19 = pkin(8) + r_i_i_C(1);
	t10 = cos(qJ(4));
	t15 = r_i_i_C(3) + qJ(5);
	t20 = pkin(4) - r_i_i_C(2);
	t7 = sin(qJ(4));
	t23 = t20 * t10 + t15 * t7 + pkin(3);
	t8 = sin(qJ(3));
	t26 = t23 * t11 + t19 * t8;
	t24 = t19 * t11;
	t9 = sin(qJ(1));
	t22 = t9 * t7;
	t21 = pkin(1) + pkin(7);
	t12 = cos(qJ(1));
	t18 = t12 * t7;
	t17 = t9 * t10;
	t16 = t12 * t10;
	t14 = pkin(3) * t8 + qJ(2) - t24;
	t4 = t8 * t16 - t22;
	t3 = t8 * t18 + t17;
	t2 = t8 * t17 + t18;
	t1 = t8 * t22 - t16;
	t5 = [t14 * t12 + t15 * t3 + t20 * t4 - t21 * t9, t9, t26 * t9, -t20 * t1 + t15 * t2, t1, 0; t15 * t1 + t21 * t12 + t14 * t9 + t20 * t2, -t12, -t26 * t12, -t15 * t4 + t20 * t3, -t3, 0; 0, 0, -t23 * t8 + t24, (t15 * t10 - t20 * t7) * t11, t11 * t7, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:48
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (83->27), mult. (186->34), div. (0->0), fcn. (218->6), ass. (0->23)
	t11 = cos(qJ(3));
	t16 = pkin(5) + pkin(8) + r_i_i_C(1);
	t10 = cos(qJ(4));
	t15 = pkin(4) + r_i_i_C(3) + qJ(6);
	t17 = r_i_i_C(2) + qJ(5);
	t7 = sin(qJ(4));
	t23 = t15 * t10 + t17 * t7 + pkin(3);
	t8 = sin(qJ(3));
	t26 = t23 * t11 + t16 * t8;
	t24 = t16 * t11;
	t9 = sin(qJ(1));
	t22 = t9 * t7;
	t21 = pkin(1) + pkin(7);
	t12 = cos(qJ(1));
	t20 = t12 * t7;
	t19 = t9 * t10;
	t18 = t12 * t10;
	t14 = pkin(3) * t8 + qJ(2) - t24;
	t4 = t8 * t18 - t22;
	t3 = t8 * t20 + t19;
	t2 = t8 * t19 + t20;
	t1 = t8 * t22 - t18;
	t5 = [t14 * t12 + t15 * t4 + t17 * t3 - t21 * t9, t9, t26 * t9, -t15 * t1 + t17 * t2, t1, t2; t17 * t1 + t21 * t12 + t14 * t9 + t15 * t2, -t12, -t26 * t12, t15 * t3 - t17 * t4, -t3, -t4; 0, 0, -t23 * t8 + t24, (t17 * t10 - t15 * t7) * t11, t11 * t7, t11 * t10;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end