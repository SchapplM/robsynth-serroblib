% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PPRRP3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRRP3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(pkin(7)), 0, 0, 0; 0, -cos(pkin(7)), 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (8->8), mult. (20->15), div. (0->0), fcn. (26->6), ass. (0->10)
	t2 = sin(pkin(7));
	t5 = sin(qJ(3));
	t10 = t2 * t5;
	t6 = cos(qJ(3));
	t9 = t2 * t6;
	t4 = cos(pkin(7));
	t8 = t4 * t5;
	t7 = t4 * t6;
	t3 = cos(pkin(8));
	t1 = [0, t2, (-t3 * t8 + t9) * r_i_i_C(1) + (-t3 * t7 - t10) * r_i_i_C(2), 0, 0; 0, -t4, (-t3 * t10 - t7) * r_i_i_C(1) + (-t3 * t9 + t8) * r_i_i_C(2), 0, 0; 1, 0, (-r_i_i_C(1) * t5 - r_i_i_C(2) * t6) * sin(pkin(8)), 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (36->20), mult. (92->38), div. (0->0), fcn. (114->8), ass. (0->20)
	t5 = sin(pkin(8));
	t9 = sin(qJ(4));
	t21 = t5 * t9;
	t20 = pkin(6) + r_i_i_C(3);
	t11 = cos(qJ(4));
	t19 = t11 * t5;
	t12 = cos(qJ(3));
	t18 = t12 * t5;
	t6 = sin(pkin(7));
	t17 = t12 * t6;
	t8 = cos(pkin(7));
	t16 = t12 * t8;
	t10 = sin(qJ(3));
	t15 = t6 * t10;
	t14 = t8 * t10;
	t13 = r_i_i_C(1) * t11 - r_i_i_C(2) * t9 + pkin(3);
	t7 = cos(pkin(8));
	t4 = t7 * t16 + t15;
	t2 = t7 * t17 - t14;
	t1 = [0, t6, t20 * t4 + t13 * (-t7 * t14 + t17), (t8 * t19 - t4 * t9) * r_i_i_C(1) + (-t11 * t4 - t8 * t21) * r_i_i_C(2), 0; 0, -t8, t20 * t2 + t13 * (-t7 * t15 - t16), (t6 * t19 - t2 * t9) * r_i_i_C(1) + (-t11 * t2 - t6 * t21) * r_i_i_C(2), 0; 1, 0, (-t13 * t10 + t20 * t12) * t5, (-t7 * t11 - t9 * t18) * r_i_i_C(1) + (-t11 * t18 + t7 * t9) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (67->22), mult. (169->38), div. (0->0), fcn. (213->8), ass. (0->25)
	t17 = sin(qJ(4));
	t19 = cos(qJ(4));
	t22 = r_i_i_C(3) + qJ(5);
	t31 = pkin(4) + r_i_i_C(1);
	t21 = t22 * t17 + t31 * t19 + pkin(3);
	t30 = pkin(6) + r_i_i_C(2);
	t13 = sin(pkin(8));
	t29 = t13 * t17;
	t28 = t13 * t19;
	t20 = cos(qJ(3));
	t27 = t13 * t20;
	t14 = sin(pkin(7));
	t18 = sin(qJ(3));
	t26 = t14 * t18;
	t25 = t14 * t20;
	t16 = cos(pkin(7));
	t24 = t16 * t18;
	t23 = t16 * t20;
	t15 = cos(pkin(8));
	t9 = t15 * t19 + t17 * t27;
	t8 = t15 * t23 + t26;
	t6 = t15 * t25 - t24;
	t3 = -t16 * t28 + t17 * t8;
	t1 = -t14 * t28 + t17 * t6;
	t2 = [0, t14, t30 * t8 + t21 * (-t15 * t24 + t25), t22 * (t16 * t29 + t19 * t8) - t31 * t3, t3; 0, -t16, t30 * t6 + t21 * (-t15 * t26 - t23), t22 * (t14 * t29 + t19 * t6) - t31 * t1, t1; 1, 0, (-t21 * t18 + t30 * t20) * t13, -t31 * t9 + t22 * (-t15 * t17 + t19 * t27), t9;];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,5);
end