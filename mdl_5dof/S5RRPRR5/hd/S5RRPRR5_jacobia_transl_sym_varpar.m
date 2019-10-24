% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:49
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRPRR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:28
	% EndTime: 2019-10-24 10:49:28
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:28
	% EndTime: 2019-10-24 10:49:28
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
	% StartTime: 2019-10-24 10:49:28
	% EndTime: 2019-10-24 10:49:28
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
	% StartTime: 2019-10-24 10:49:28
	% EndTime: 2019-10-24 10:49:28
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (40->8), mult. (30->8), div. (0->0), fcn. (32->6), ass. (0->8)
	t15 = r_i_i_C(3) + qJ(3);
	t14 = r_i_i_C(2) * sin(pkin(9)) - r_i_i_C(1) * cos(pkin(9)) - pkin(2);
	t7 = qJ(1) + qJ(2);
	t5 = sin(t7);
	t6 = cos(t7);
	t11 = t14 * t5 + t15 * t6;
	t10 = t14 * t6 - t15 * t5;
	t1 = [0, 0, 0, 0, 0; -cos(qJ(1)) * pkin(1) + t10, t10, t6, 0, 0; -sin(qJ(1)) * pkin(1) + t11, t11, t5, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:28
	% EndTime: 2019-10-24 10:49:28
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (69->12), mult. (44->13), div. (0->0), fcn. (46->7), ass. (0->13)
	t10 = pkin(9) + qJ(4);
	t6 = sin(t10);
	t7 = cos(t10);
	t21 = t7 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t20 = r_i_i_C(3) + pkin(7) + qJ(3);
	t19 = -cos(pkin(9)) * pkin(3) - pkin(2) - t21;
	t15 = r_i_i_C(1) * t6 + r_i_i_C(2) * t7;
	t11 = qJ(1) + qJ(2);
	t8 = sin(t11);
	t9 = cos(t11);
	t14 = t19 * t9 - t20 * t8;
	t13 = t19 * t8 + t20 * t9;
	t1 = [0, 0, 0, t21, 0; -cos(qJ(1)) * pkin(1) + t14, t14, t9, t15 * t8, 0; -sin(qJ(1)) * pkin(1) + t13, t13, t8, -t15 * t9, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:28
	% EndTime: 2019-10-24 10:49:28
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (122->19), mult. (63->20), div. (0->0), fcn. (65->9), ass. (0->21)
	t29 = r_i_i_C(3) + pkin(8) + pkin(7) + qJ(3);
	t17 = pkin(9) + qJ(4);
	t13 = qJ(5) + t17;
	t9 = sin(t13);
	t28 = t9 * r_i_i_C(2);
	t18 = qJ(1) + qJ(2);
	t14 = sin(t18);
	t24 = t14 * t9;
	t10 = cos(t13);
	t25 = r_i_i_C(2) * t10;
	t27 = r_i_i_C(1) * t24 + t14 * t25;
	t26 = pkin(4) * sin(t17);
	t7 = t10 * r_i_i_C(1);
	t23 = t7 - t28;
	t8 = pkin(4) * cos(t17);
	t22 = -t8 - cos(pkin(9)) * pkin(3) - pkin(2) - t7;
	t21 = -r_i_i_C(1) * t9 - t25;
	t15 = cos(t18);
	t20 = (t22 + t28) * t15 - t29 * t14;
	t19 = r_i_i_C(2) * t24 + t22 * t14 + t29 * t15;
	t1 = [0, 0, 0, t23 + t8, t23; -cos(qJ(1)) * pkin(1) + t20, t20, t15, t14 * t26 + t27, t27; -sin(qJ(1)) * pkin(1) + t19, t19, t14, (t21 - t26) * t15, t21 * t15;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end