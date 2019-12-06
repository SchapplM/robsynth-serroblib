% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRRR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PPRRR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRRR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:15:17
	% EndTime: 2019-12-05 15:15:17
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:15:17
	% EndTime: 2019-12-05 15:15:17
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:15:17
	% EndTime: 2019-12-05 15:15:17
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(pkin(8)), 0, 0, 0; 0, -cos(pkin(8)), 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:15:17
	% EndTime: 2019-12-05 15:15:17
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (10->4), mult. (10->6), div. (0->0), fcn. (12->4), ass. (0->7)
	t3 = pkin(9) + qJ(3);
	t1 = sin(t3);
	t2 = cos(t3);
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t2;
	t5 = cos(pkin(8));
	t4 = sin(pkin(8));
	t7 = [0, t4, t6 * t5, 0, 0; 0, -t5, t6 * t4, 0, 0; 1, 0, t2 * r_i_i_C(1) - t1 * r_i_i_C(2), 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:15:17
	% EndTime: 2019-12-05 15:15:17
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (41->14), mult. (51->23), div. (0->0), fcn. (57->6), ass. (0->15)
	t4 = sin(pkin(8));
	t6 = sin(qJ(4));
	t14 = t4 * t6;
	t7 = cos(qJ(4));
	t13 = t4 * t7;
	t5 = cos(pkin(8));
	t12 = t5 * t6;
	t11 = t5 * t7;
	t10 = pkin(6) + r_i_i_C(3);
	t9 = t7 * r_i_i_C(1) - t6 * r_i_i_C(2) + pkin(3);
	t3 = pkin(9) + qJ(3);
	t1 = sin(t3);
	t2 = cos(t3);
	t8 = -t9 * t1 + t10 * t2;
	t15 = [0, t4, t8 * t5, (-t2 * t12 + t13) * r_i_i_C(1) + (-t2 * t11 - t14) * r_i_i_C(2), 0; 0, -t5, t8 * t4, (-t2 * t14 - t11) * r_i_i_C(1) + (-t2 * t13 + t12) * r_i_i_C(2), 0; 1, 0, t10 * t1 + t9 * t2, (-r_i_i_C(1) * t6 - r_i_i_C(2) * t7) * t1, 0;];
	Ja_transl = t15;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:15:17
	% EndTime: 2019-12-05 15:15:17
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (94->21), mult. (86->33), div. (0->0), fcn. (96->8), ass. (0->22)
	t13 = cos(pkin(8));
	t11 = qJ(4) + qJ(5);
	t9 = cos(t11);
	t22 = t13 * t9;
	t8 = sin(t11);
	t23 = t13 * t8;
	t12 = sin(pkin(8));
	t24 = t12 * t9;
	t25 = t12 * t8;
	t10 = pkin(9) + qJ(3);
	t7 = cos(t10);
	t27 = (-t7 * t25 - t22) * r_i_i_C(1) + (-t7 * t24 + t23) * r_i_i_C(2);
	t26 = (-t7 * t23 + t24) * r_i_i_C(1) + (-t7 * t22 - t25) * r_i_i_C(2);
	t14 = sin(qJ(4));
	t21 = t14 * t7;
	t20 = r_i_i_C(3) + pkin(7) + pkin(6);
	t19 = -r_i_i_C(1) * t8 - r_i_i_C(2) * t9;
	t15 = cos(qJ(4));
	t18 = pkin(4) * t15 + r_i_i_C(1) * t9 - r_i_i_C(2) * t8 + pkin(3);
	t6 = sin(t10);
	t17 = -t18 * t6 + t20 * t7;
	t1 = [0, t12, t17 * t13, (t12 * t15 - t13 * t21) * pkin(4) + t26, t26; 0, -t13, t17 * t12, (-t12 * t21 - t13 * t15) * pkin(4) + t27, t27; 1, 0, t18 * t7 + t20 * t6, (-pkin(4) * t14 + t19) * t6, t19 * t6;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end