% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRPRP6
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRPRP6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRPRP6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_jacobia_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:41:37
	% EndTime: 2019-12-05 15:41:37
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:41:37
	% EndTime: 2019-12-05 15:41:37
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:41:37
	% EndTime: 2019-12-05 15:41:38
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->2), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t3 = sin(qJ(2));
	t4 = cos(qJ(2));
	t5 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t1 = [0, t5 * cos(pkin(7)), 0, 0, 0; 0, t5 * sin(pkin(7)), 0, 0, 0; 1, t4 * r_i_i_C(1) - t3 * r_i_i_C(2), 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:41:37
	% EndTime: 2019-12-05 15:41:38
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (10->5), mult. (22->8), div. (0->0), fcn. (25->4), ass. (0->8)
	t7 = pkin(2) - r_i_i_C(2);
	t6 = r_i_i_C(3) + qJ(3);
	t3 = sin(qJ(2));
	t4 = cos(qJ(2));
	t5 = -t3 * t7 + t4 * t6;
	t2 = cos(pkin(7));
	t1 = sin(pkin(7));
	t8 = [0, t5 * t2, t2 * t3, 0, 0; 0, t5 * t1, t1 * t3, 0, 0; 1, t3 * t6 + t4 * t7, -t4, 0, 0;];
	Ja_transl = t8;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:41:37
	% EndTime: 2019-12-05 15:41:38
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (23->13), mult. (58->27), div. (0->0), fcn. (65->6), ass. (0->12)
	t3 = sin(qJ(4));
	t4 = sin(qJ(2));
	t11 = t3 * t4;
	t5 = cos(qJ(4));
	t10 = t4 * t5;
	t9 = pkin(2) + pkin(6) + r_i_i_C(3);
	t8 = t3 * r_i_i_C(1) + t5 * r_i_i_C(2) + qJ(3);
	t6 = cos(qJ(2));
	t7 = -t9 * t4 + t8 * t6;
	t2 = cos(pkin(7));
	t1 = sin(pkin(7));
	t12 = [0, t7 * t2, t2 * t4, (-t1 * t3 + t2 * t10) * r_i_i_C(1) + (-t1 * t5 - t2 * t11) * r_i_i_C(2), 0; 0, t7 * t1, t1 * t4, (t1 * t10 + t2 * t3) * r_i_i_C(1) + (-t1 * t11 + t2 * t5) * r_i_i_C(2), 0; 1, t8 * t4 + t9 * t6, -t6, (-r_i_i_C(1) * t5 + r_i_i_C(2) * t3) * t6, 0;];
	Ja_transl = t12;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:41:37
	% EndTime: 2019-12-05 15:41:38
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (41->16), mult. (101->28), div. (0->0), fcn. (117->6), ass. (0->16)
	t8 = sin(qJ(4));
	t9 = sin(qJ(2));
	t18 = t8 * t9;
	t17 = pkin(4) + r_i_i_C(1);
	t10 = cos(qJ(4));
	t16 = t10 * t9;
	t15 = r_i_i_C(3) + qJ(5);
	t14 = pkin(2) + pkin(6) + r_i_i_C(2);
	t13 = -t15 * t10 + t17 * t8 + qJ(3);
	t11 = cos(qJ(2));
	t12 = t13 * t11 - t14 * t9;
	t7 = cos(pkin(7));
	t6 = sin(pkin(7));
	t3 = t6 * t16 + t7 * t8;
	t1 = -t7 * t16 + t6 * t8;
	t2 = [0, t12 * t7, t7 * t9, t15 * (t10 * t6 + t7 * t18) - t17 * t1, t1; 0, t12 * t6, t6 * t9, -t15 * (t10 * t7 - t6 * t18) + t17 * t3, -t3; 1, t14 * t11 + t13 * t9, -t11, (-t17 * t10 - t15 * t8) * t11, t11 * t10;];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,5);
end