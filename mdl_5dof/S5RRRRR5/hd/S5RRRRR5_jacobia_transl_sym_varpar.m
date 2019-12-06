% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRRR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:59:11
	% EndTime: 2019-12-05 18:59:11
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:59:11
	% EndTime: 2019-12-05 18:59:11
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [0, 0, 0, 0, 0; -r_i_i_C(1) * t2 + r_i_i_C(2) * t1, 0, 0, 0, 0; -r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:59:11
	% EndTime: 2019-12-05 18:59:11
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (14->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->6)
	t4 = qJ(1) + qJ(2);
	t2 = sin(t4);
	t3 = cos(t4);
	t6 = -t3 * r_i_i_C(1) + t2 * r_i_i_C(2);
	t5 = -t2 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t1 = [0, 0, 0, 0, 0; -cos(qJ(1)) * pkin(1) + t6, t6, 0, 0, 0; -sin(qJ(1)) * pkin(1) + t5, t5, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:59:11
	% EndTime: 2019-12-05 18:59:11
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (40->8), mult. (18->8), div. (0->0), fcn. (18->6), ass. (0->9)
	t7 = qJ(1) + qJ(2);
	t6 = qJ(3) + t7;
	t2 = sin(t6);
	t3 = cos(t6);
	t11 = -t3 * r_i_i_C(1) + t2 * r_i_i_C(2);
	t10 = -t2 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t9 = -pkin(2) * cos(t7) + t11;
	t8 = -pkin(2) * sin(t7) + t10;
	t1 = [0, 0, 0, 0, 0; -cos(qJ(1)) * pkin(1) + t9, t9, t11, 0, 0; -sin(qJ(1)) * pkin(1) + t8, t8, t10, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:59:11
	% EndTime: 2019-12-05 18:59:11
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (105->13), mult. (58->14), div. (0->0), fcn. (58->8), ass. (0->15)
	t11 = sin(qJ(4));
	t12 = cos(qJ(4));
	t23 = t12 * r_i_i_C(1) - t11 * r_i_i_C(2);
	t22 = pkin(8) + r_i_i_C(3);
	t21 = -pkin(3) - t23;
	t10 = qJ(1) + qJ(2);
	t17 = r_i_i_C(1) * t11 + r_i_i_C(2) * t12;
	t9 = qJ(3) + t10;
	t5 = sin(t9);
	t6 = cos(t9);
	t16 = t21 * t5 + t22 * t6;
	t15 = t21 * t6 - t22 * t5;
	t14 = -pkin(2) * sin(t10) + t16;
	t13 = -pkin(2) * cos(t10) + t15;
	t1 = [0, 0, 0, t23, 0; -cos(qJ(1)) * pkin(1) + t13, t13, t15, t17 * t5, 0; -sin(qJ(1)) * pkin(1) + t14, t14, t16, -t17 * t6, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:59:11
	% EndTime: 2019-12-05 18:59:11
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (159->19), mult. (79->21), div. (0->0), fcn. (79->10), ass. (0->23)
	t33 = r_i_i_C(3) + pkin(9) + pkin(8);
	t17 = qJ(4) + qJ(5);
	t11 = sin(t17);
	t18 = qJ(1) + qJ(2);
	t15 = qJ(3) + t18;
	t8 = sin(t15);
	t28 = t11 * t8;
	t13 = cos(t17);
	t29 = r_i_i_C(2) * t13;
	t32 = r_i_i_C(1) * t28 + t8 * t29;
	t31 = sin(qJ(4)) * pkin(4);
	t30 = r_i_i_C(2) * t11;
	t7 = t13 * r_i_i_C(1);
	t27 = t7 - t30;
	t16 = cos(qJ(4)) * pkin(4);
	t26 = -t16 - pkin(3) - t7;
	t25 = -r_i_i_C(1) * t11 - t29;
	t9 = cos(t15);
	t24 = (t26 + t30) * t9 - t33 * t8;
	t23 = r_i_i_C(2) * t28 + t26 * t8 + t33 * t9;
	t22 = -pkin(2) * cos(t18) + t24;
	t21 = -pkin(2) * sin(t18) + t23;
	t1 = [0, 0, 0, t16 + t27, t27; -cos(qJ(1)) * pkin(1) + t22, t22, t24, t8 * t31 + t32, t32; -sin(qJ(1)) * pkin(1) + t21, t21, t23, (t25 - t31) * t9, t25 * t9;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end