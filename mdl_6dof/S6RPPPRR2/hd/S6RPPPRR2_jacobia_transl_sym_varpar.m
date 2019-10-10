% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:27
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPPRR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
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
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
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
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (19->8), mult. (10->6), div. (0->0), fcn. (12->4), ass. (0->6)
	t5 = pkin(2) - r_i_i_C(2);
	t4 = r_i_i_C(3) + qJ(3);
	t3 = qJ(1) + pkin(9);
	t2 = cos(t3);
	t1 = sin(t3);
	t6 = [-sin(qJ(1)) * pkin(1) + t4 * t2 - t5 * t1, 0, t1, 0, 0, 0; cos(qJ(1)) * pkin(1) + t5 * t2 + t4 * t1, 0, -t2, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_transl = t6;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (29->9), mult. (18->8), div. (0->0), fcn. (22->6), ass. (0->6)
	t7 = pkin(2) + r_i_i_C(3) + qJ(4);
	t6 = r_i_i_C(1) * sin(pkin(10)) + r_i_i_C(2) * cos(pkin(10)) + qJ(3);
	t3 = qJ(1) + pkin(9);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [-sin(qJ(1)) * pkin(1) + t6 * t2 - t7 * t1, 0, t1, t2, 0, 0; cos(qJ(1)) * pkin(1) + t7 * t2 + t6 * t1, 0, -t2, t1, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (52->13), mult. (32->13), div. (0->0), fcn. (36->7), ass. (0->11)
	t12 = pkin(2) + r_i_i_C(3) + pkin(7) + qJ(4);
	t5 = pkin(10) + qJ(5);
	t1 = sin(t5);
	t3 = cos(t5);
	t11 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t10 = -t1 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t9 = pkin(4) * sin(pkin(10)) + qJ(3) - t10;
	t6 = qJ(1) + pkin(9);
	t4 = cos(t6);
	t2 = sin(t6);
	t7 = [-sin(qJ(1)) * pkin(1) - t12 * t2 + t9 * t4, 0, t2, t4, t11 * t2, 0; cos(qJ(1)) * pkin(1) + t12 * t4 + t9 * t2, 0, -t4, t2, -t11 * t4, 0; 0, 1, 0, 0, t10, 0;];
	Ja_transl = t7;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (119->28), mult. (91->35), div. (0->0), fcn. (103->9), ass. (0->23)
	t22 = pkin(8) + r_i_i_C(3);
	t9 = pkin(10) + qJ(5);
	t7 = cos(t9);
	t24 = t22 * t7;
	t13 = sin(qJ(6));
	t14 = cos(qJ(6));
	t16 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + pkin(5);
	t5 = sin(t9);
	t23 = t16 * t7 + t22 * t5;
	t10 = qJ(1) + pkin(9);
	t6 = sin(t10);
	t21 = t14 * t6;
	t8 = cos(t10);
	t20 = t14 * t8;
	t19 = t6 * t13;
	t18 = t8 * t13;
	t17 = pkin(2) + pkin(7) + qJ(4);
	t15 = pkin(4) * sin(pkin(10)) + t5 * pkin(5) - t24 + qJ(3);
	t4 = t5 * t20 - t19;
	t3 = t5 * t18 + t21;
	t2 = t5 * t21 + t18;
	t1 = -t5 * t19 + t20;
	t11 = [-sin(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t17 * t6 + t15 * t8, 0, t6, t8, t23 * t6, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; cos(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t17 * t8 + t15 * t6, 0, -t8, t6, -t23 * t8, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 1, 0, 0, -t16 * t5 + t24, (-r_i_i_C(1) * t13 - r_i_i_C(2) * t14) * t7;];
	Ja_transl = t11;
else
	Ja_transl=NaN(3,6);
end