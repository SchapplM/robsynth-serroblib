% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRR6
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
% Datum: 2019-10-24 10:52
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRRR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:49
	% EndTime: 2019-10-24 10:52:49
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:49
	% EndTime: 2019-10-24 10:52:49
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
	% StartTime: 2019-10-24 10:52:49
	% EndTime: 2019-10-24 10:52:49
	% DurationCPUTime: 0.06s
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
	% StartTime: 2019-10-24 10:52:49
	% EndTime: 2019-10-24 10:52:49
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-24 10:52:49
	% EndTime: 2019-10-24 10:52:49
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (85->16), mult. (59->19), div. (0->0), fcn. (59->8), ass. (0->20)
	t28 = r_i_i_C(3) + pkin(8) + pkin(7);
	t14 = qJ(3) + qJ(4);
	t9 = sin(t14);
	t27 = t9 * r_i_i_C(2);
	t15 = qJ(1) + qJ(2);
	t10 = sin(t15);
	t23 = t10 * t9;
	t11 = cos(t14);
	t24 = r_i_i_C(2) * t11;
	t26 = r_i_i_C(1) * t23 + t10 * t24;
	t25 = sin(qJ(3)) * pkin(3);
	t6 = t11 * r_i_i_C(1);
	t22 = t6 - t27;
	t13 = cos(qJ(3)) * pkin(3);
	t21 = -t13 - pkin(2) - t6;
	t20 = -r_i_i_C(1) * t9 - t24;
	t12 = cos(t15);
	t19 = (t21 + t27) * t12 - t28 * t10;
	t18 = r_i_i_C(2) * t23 + t21 * t10 + t28 * t12;
	t1 = [0, 0, t13 + t22, t22, 0; -cos(qJ(1)) * pkin(1) + t19, t19, t10 * t25 + t26, t26, 0; -sin(qJ(1)) * pkin(1) + t18, t18, (t20 - t25) * t12, t20 * t12, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:49
	% EndTime: 2019-10-24 10:52:49
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (150->21), mult. (81->23), div. (0->0), fcn. (81->10), ass. (0->24)
	t32 = r_i_i_C(3) + pkin(9) + pkin(8) + pkin(7);
	t20 = qJ(1) + qJ(2);
	t14 = sin(t20);
	t19 = qJ(3) + qJ(4);
	t16 = qJ(5) + t19;
	t11 = sin(t16);
	t27 = t11 * t14;
	t12 = cos(t16);
	t28 = r_i_i_C(2) * t12;
	t31 = r_i_i_C(1) * t27 + t14 * t28;
	t30 = pkin(4) * sin(t19);
	t29 = r_i_i_C(2) * t11;
	t8 = t12 * r_i_i_C(1);
	t10 = pkin(4) * cos(t19);
	t17 = cos(qJ(3)) * pkin(3);
	t26 = -t10 - t17 - pkin(2) - t8;
	t25 = t8 - t29;
	t24 = t10 + t25;
	t23 = -r_i_i_C(1) * t11 - t28;
	t15 = cos(t20);
	t22 = (t26 + t29) * t15 - t32 * t14;
	t21 = r_i_i_C(2) * t27 + t26 * t14 + t32 * t15;
	t6 = -t30 - sin(qJ(3)) * pkin(3);
	t1 = [0, 0, t17 + t24, t24, t25; -cos(qJ(1)) * pkin(1) + t22, t22, -t14 * t6 + t31, t14 * t30 + t31, t31; -sin(qJ(1)) * pkin(1) + t21, t21, (t23 + t6) * t15, (t23 - t30) * t15, t23 * t15;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end