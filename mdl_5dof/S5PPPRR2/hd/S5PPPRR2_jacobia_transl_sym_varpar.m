% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPPRR2
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
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PPPRR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPPRR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPPRR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:52
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:52
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:52
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(pkin(7)), 0, 0, 0; 0, -cos(pkin(7)), 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:52
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (2->2), div. (0->0), fcn. (7->4), ass. (0->4)
	t3 = cos(pkin(7));
	t2 = sin(pkin(7));
	t1 = sin(pkin(8));
	t4 = [0, t2, t3 * t1, 0, 0; 0, -t3, t2 * t1, 0, 0; 1, 0, -cos(pkin(8)), 0, 0;];
	Ja_transl = t4;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:52
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (15->13), mult. (38->27), div. (0->0), fcn. (53->8), ass. (0->15)
	t4 = sin(pkin(8));
	t5 = sin(pkin(7));
	t14 = t5 * t4;
	t6 = cos(pkin(9));
	t7 = cos(pkin(8));
	t13 = t6 * t7;
	t8 = cos(pkin(7));
	t12 = t8 * t4;
	t10 = cos(qJ(4));
	t11 = t10 * t4;
	t9 = sin(qJ(4));
	t3 = sin(pkin(9));
	t2 = t8 * t13 + t5 * t3;
	t1 = t5 * t13 - t8 * t3;
	t15 = [0, t5, t12, (t8 * t11 - t2 * t9) * r_i_i_C(1) + (-t2 * t10 - t9 * t12) * r_i_i_C(2), 0; 0, -t8, t14, (-t1 * t9 + t11 * t5) * r_i_i_C(1) + (-t1 * t10 - t9 * t14) * r_i_i_C(2), 0; 1, 0, -t7, (-t4 * t6 * t9 - t7 * t10) * r_i_i_C(1) + (-t11 * t6 + t7 * t9) * r_i_i_C(2), 0;];
	Ja_transl = t15;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:52
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (62->27), mult. (160->52), div. (0->0), fcn. (210->10), ass. (0->26)
	t28 = pkin(6) + r_i_i_C(3);
	t11 = sin(pkin(9));
	t12 = sin(pkin(8));
	t27 = t11 * t12;
	t14 = cos(pkin(9));
	t26 = t12 * t14;
	t13 = sin(pkin(7));
	t25 = t13 * t12;
	t15 = cos(pkin(8));
	t24 = t13 * t15;
	t16 = cos(pkin(7));
	t23 = t15 * t16;
	t22 = t16 * t12;
	t17 = sin(qJ(5));
	t19 = cos(qJ(5));
	t21 = r_i_i_C(1) * t19 - r_i_i_C(2) * t17 + pkin(4);
	t20 = cos(qJ(4));
	t18 = sin(qJ(4));
	t10 = -t15 * t18 + t20 * t26;
	t8 = t11 * t13 + t14 * t23;
	t7 = t11 * t23 - t13 * t14;
	t6 = -t11 * t16 + t14 * t24;
	t5 = t11 * t24 + t14 * t16;
	t4 = t18 * t22 + t20 * t8;
	t2 = t18 * t25 + t20 * t6;
	t1 = [0, t13, t22, t28 * t4 + t21 * (-t8 * t18 + t20 * t22), (-t17 * t4 + t19 * t7) * r_i_i_C(1) + (-t17 * t7 - t19 * t4) * r_i_i_C(2); 0, -t16, t25, t28 * t2 + t21 * (-t6 * t18 + t20 * t25), (-t17 * t2 + t19 * t5) * r_i_i_C(1) + (-t17 * t5 - t19 * t2) * r_i_i_C(2); 1, 0, -t15, t28 * t10 + t21 * (-t15 * t20 - t18 * t26), (-t10 * t17 + t19 * t27) * r_i_i_C(1) + (-t10 * t19 - t17 * t27) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end