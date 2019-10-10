% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:34
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRPR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
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
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (23->9), mult. (16->8), div. (0->0), fcn. (18->6), ass. (0->6)
	t7 = r_i_i_C(3) + qJ(3);
	t6 = r_i_i_C(1) * cos(pkin(10)) - r_i_i_C(2) * sin(pkin(10)) + pkin(2);
	t3 = qJ(1) + pkin(9);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [-sin(qJ(1)) * pkin(1) + t7 * t2 - t6 * t1, 0, t1, 0, 0, 0; cos(qJ(1)) * pkin(1) + t7 * t1 + t6 * t2, 0, -t2, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (44->12), mult. (28->13), div. (0->0), fcn. (30->7), ass. (0->11)
	t12 = r_i_i_C(3) + pkin(7) + qJ(3);
	t6 = pkin(10) + qJ(4);
	t2 = sin(t6);
	t4 = cos(t6);
	t11 = r_i_i_C(1) * t4 - r_i_i_C(2) * t2;
	t10 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t4;
	t9 = cos(pkin(10)) * pkin(3) + pkin(2) + t11;
	t7 = qJ(1) + pkin(9);
	t5 = cos(t7);
	t3 = sin(t7);
	t1 = [-sin(qJ(1)) * pkin(1) + t12 * t5 - t9 * t3, 0, t3, t10 * t5, 0, 0; cos(qJ(1)) * pkin(1) + t12 * t3 + t9 * t5, 0, -t5, t10 * t3, 0, 0; 0, 1, 0, t11, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (96->18), mult. (69->19), div. (0->0), fcn. (78->9), ass. (0->15)
	t8 = sin(pkin(11));
	t9 = cos(pkin(11));
	t14 = r_i_i_C(1) * t9 - r_i_i_C(2) * t8 + pkin(4);
	t15 = r_i_i_C(3) + qJ(5);
	t6 = pkin(10) + qJ(4);
	t2 = sin(t6);
	t4 = cos(t6);
	t12 = t14 * t4 + t15 * t2;
	t16 = cos(pkin(10)) * pkin(3) + pkin(2) + t12;
	t13 = r_i_i_C(1) * t8 + r_i_i_C(2) * t9 + pkin(7) + qJ(3);
	t11 = -t14 * t2 + t15 * t4;
	t7 = qJ(1) + pkin(9);
	t5 = cos(t7);
	t3 = sin(t7);
	t1 = [-sin(qJ(1)) * pkin(1) + t13 * t5 - t16 * t3, 0, t3, t11 * t5, t5 * t2, 0; cos(qJ(1)) * pkin(1) + t13 * t3 + t16 * t5, 0, -t5, t11 * t3, t3 * t2, 0; 0, 1, 0, t12, -t4, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (155->31), mult. (98->41), div. (0->0), fcn. (111->11), ass. (0->23)
	t14 = pkin(10) + qJ(4);
	t11 = cos(t14);
	t25 = r_i_i_C(3) + pkin(8) + qJ(5);
	t8 = sin(t14);
	t23 = t25 * t8;
	t5 = cos(pkin(11)) * pkin(5) + pkin(4);
	t27 = t23 + t11 * t5 + cos(pkin(10)) * pkin(3) + pkin(2);
	t15 = qJ(1) + pkin(9);
	t9 = sin(t15);
	t26 = t11 * t9;
	t12 = cos(t15);
	t24 = t11 * t12;
	t21 = pkin(5) * sin(pkin(11)) + pkin(7) + qJ(3);
	t13 = pkin(11) + qJ(6);
	t10 = cos(t13);
	t7 = sin(t13);
	t20 = r_i_i_C(1) * t10 - r_i_i_C(2) * t7 + t5;
	t19 = t25 * t11 - t20 * t8;
	t4 = t10 * t24 + t7 * t9;
	t3 = t10 * t9 - t7 * t24;
	t2 = -t10 * t26 + t12 * t7;
	t1 = t10 * t12 + t7 * t26;
	t6 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t21 * t12 - t27 * t9, 0, t9, t19 * t12, t12 * t8, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t21 * t9 + t27 * t12, 0, -t12, t19 * t9, t9 * t8, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, 1, 0, t20 * t11 + t23, -t11, (-r_i_i_C(1) * t7 - r_i_i_C(2) * t10) * t8;];
	Ja_transl = t6;
else
	Ja_transl=NaN(3,6);
end