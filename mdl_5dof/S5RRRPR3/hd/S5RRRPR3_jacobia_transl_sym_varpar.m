% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRPR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:43:51
	% EndTime: 2019-12-05 18:43:51
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:43:51
	% EndTime: 2019-12-05 18:43:51
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
	% StartTime: 2019-12-05 18:43:51
	% EndTime: 2019-12-05 18:43:51
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-12-05 18:43:51
	% EndTime: 2019-12-05 18:43:51
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (45->10), mult. (40->12), div. (0->0), fcn. (40->6), ass. (0->12)
	t8 = sin(qJ(3));
	t9 = cos(qJ(3));
	t18 = t9 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t17 = pkin(7) + r_i_i_C(3);
	t16 = -pkin(2) - t18;
	t12 = r_i_i_C(1) * t8 + r_i_i_C(2) * t9;
	t7 = qJ(1) + qJ(2);
	t5 = sin(t7);
	t6 = cos(t7);
	t11 = t16 * t5 + t17 * t6;
	t10 = t16 * t6 - t17 * t5;
	t1 = [0, 0, t18, 0, 0; -cos(qJ(1)) * pkin(1) + t10, t10, t12 * t5, 0, 0; -sin(qJ(1)) * pkin(1) + t11, t11, -t12 * t6, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:43:51
	% EndTime: 2019-12-05 18:43:51
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (74->13), mult. (49->14), div. (0->0), fcn. (51->8), ass. (0->13)
	t11 = qJ(3) + pkin(9);
	t6 = sin(t11);
	t7 = cos(t11);
	t24 = -t6 * r_i_i_C(2) + cos(qJ(3)) * pkin(3) + t7 * r_i_i_C(1);
	t23 = r_i_i_C(3) + qJ(4) + pkin(7);
	t22 = -pkin(2) - t24;
	t17 = sin(qJ(3)) * pkin(3) + r_i_i_C(1) * t6 + r_i_i_C(2) * t7;
	t12 = qJ(1) + qJ(2);
	t8 = sin(t12);
	t9 = cos(t12);
	t16 = t22 * t9 - t23 * t8;
	t15 = t22 * t8 + t23 * t9;
	t1 = [0, 0, t24, 0, 0; -cos(qJ(1)) * pkin(1) + t16, t16, t17 * t8, t9, 0; -sin(qJ(1)) * pkin(1) + t15, t15, -t17 * t9, t8, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:43:51
	% EndTime: 2019-12-05 18:43:51
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (125->19), mult. (66->21), div. (0->0), fcn. (68->10), ass. (0->21)
	t30 = r_i_i_C(3) + pkin(8) + qJ(4) + pkin(7);
	t19 = qJ(1) + qJ(2);
	t14 = sin(t19);
	t18 = qJ(3) + pkin(9);
	t13 = qJ(5) + t18;
	t10 = sin(t13);
	t25 = t10 * t14;
	t11 = cos(t13);
	t28 = r_i_i_C(2) * t11;
	t29 = r_i_i_C(1) * t25 + t14 * t28;
	t27 = t10 * r_i_i_C(2);
	t8 = t11 * r_i_i_C(1);
	t26 = cos(qJ(3)) * pkin(3) + pkin(4) * cos(t18);
	t24 = -pkin(2) - t26 - t8;
	t23 = t8 - t27;
	t22 = -r_i_i_C(1) * t10 - t28;
	t15 = cos(t19);
	t21 = (t24 + t27) * t15 - t30 * t14;
	t20 = r_i_i_C(2) * t25 + t24 * t14 + t30 * t15;
	t6 = -pkin(4) * sin(t18) - sin(qJ(3)) * pkin(3);
	t1 = [0, 0, t23 + t26, 0, t23; -cos(qJ(1)) * pkin(1) + t21, t21, -t14 * t6 + t29, t15, t29; -sin(qJ(1)) * pkin(1) + t20, t20, (t22 + t6) * t15, t14, t22 * t15;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end