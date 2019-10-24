% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:41
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPPRR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPRR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:41:09
	% EndTime: 2019-10-24 10:41:09
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:41:08
	% EndTime: 2019-10-24 10:41:09
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
	% StartTime: 2019-10-24 10:41:09
	% EndTime: 2019-10-24 10:41:09
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (8->5), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = -r_i_i_C(1) * cos(pkin(8)) + r_i_i_C(2) * sin(pkin(8)) - pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t3 * t6 + t4 * t5, t4, 0, 0, 0; t3 * t5 + t4 * t6, t3, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:41:09
	% EndTime: 2019-10-24 10:41:09
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (18->11), mult. (34->12), div. (0->0), fcn. (43->6), ass. (0->9)
	t1 = sin(pkin(9));
	t3 = cos(pkin(9));
	t8 = t1 * r_i_i_C(1) + t3 * r_i_i_C(2) + qJ(2);
	t2 = sin(pkin(8));
	t4 = cos(pkin(8));
	t7 = -pkin(1) + (-r_i_i_C(3) - qJ(3)) * t2 + (-r_i_i_C(1) * t3 + r_i_i_C(2) * t1 - pkin(2)) * t4;
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t9 = [0, 0, -t4, 0, 0; -t5 * t8 + t7 * t6, t6, -t5 * t2, 0, 0; t7 * t5 + t6 * t8, t5, t6 * t2, 0, 0;];
	Ja_transl = t9;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:41:09
	% EndTime: 2019-10-24 10:41:09
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (49->22), mult. (60->31), div. (0->0), fcn. (73->8), ass. (0->16)
	t11 = cos(pkin(8));
	t13 = sin(qJ(1));
	t18 = t11 * t13;
	t14 = cos(qJ(1));
	t17 = t11 * t14;
	t16 = sin(pkin(9)) * pkin(3) + qJ(2);
	t10 = sin(pkin(8));
	t15 = -t11 * (cos(pkin(9)) * pkin(3) + pkin(2)) - pkin(1) + (-r_i_i_C(3) - pkin(6) - qJ(3)) * t10;
	t8 = pkin(9) + qJ(4);
	t7 = cos(t8);
	t6 = sin(t8);
	t4 = -t13 * t6 - t7 * t17;
	t3 = -t13 * t7 + t6 * t17;
	t2 = -t14 * t6 + t7 * t18;
	t1 = t14 * t7 + t6 * t18;
	t5 = [0, 0, -t11, (-r_i_i_C(1) * t6 - r_i_i_C(2) * t7) * t10, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t16 * t13 + t15 * t14, t14, -t13 * t10, t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0; -t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * t13 + t16 * t14, t13, t14 * t10, -t3 * r_i_i_C(1) + t4 * r_i_i_C(2), 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:41:08
	% EndTime: 2019-10-24 10:41:09
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (114->31), mult. (96->40), div. (0->0), fcn. (113->10), ass. (0->24)
	t17 = pkin(9) + qJ(4);
	t15 = qJ(5) + t17;
	t12 = cos(t15);
	t19 = cos(pkin(8));
	t21 = cos(qJ(1));
	t11 = sin(t15);
	t20 = sin(qJ(1));
	t26 = t20 * t11;
	t5 = t12 * t21 + t19 * t26;
	t25 = t20 * t12;
	t6 = -t11 * t21 + t19 * t25;
	t30 = t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t27 = t19 * t21;
	t7 = t11 * t27 - t25;
	t8 = -t12 * t27 - t26;
	t29 = -t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
	t13 = sin(t17);
	t28 = pkin(4) * t13;
	t24 = qJ(2) + t28 + sin(pkin(9)) * pkin(3);
	t23 = -r_i_i_C(1) * t11 - r_i_i_C(2) * t12;
	t14 = cos(t17);
	t18 = sin(pkin(8));
	t22 = -t19 * (pkin(4) * t14 + cos(pkin(9)) * pkin(3) + pkin(2)) - pkin(1) + (-r_i_i_C(3) - pkin(7) - pkin(6) - qJ(3)) * t18;
	t1 = [0, 0, -t19, (t23 - t28) * t18, t23 * t18; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) - t24 * t20 + t22 * t21, t21, -t20 * t18, (t13 * t19 * t20 + t14 * t21) * pkin(4) + t30, t30; -t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t22 * t20 + t24 * t21, t20, t21 * t18, (-t13 * t27 + t14 * t20) * pkin(4) + t29, t29;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end