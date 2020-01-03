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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (14->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->6)
	t6 = qJ(1) + qJ(2);
	t4 = sin(t6);
	t5 = cos(t6);
	t8 = t4 * r_i_i_C(1) + t5 * r_i_i_C(2);
	t7 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t1 = [0, 0, 0, 0, 0; cos(qJ(1)) * pkin(1) + t7, t7, 0, 0, 0; sin(qJ(1)) * pkin(1) + t8, t8, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (42->10), mult. (30->8), div. (0->0), fcn. (32->6), ass. (0->8)
	t17 = r_i_i_C(3) + qJ(3);
	t16 = -r_i_i_C(2) * sin(pkin(9)) + r_i_i_C(1) * cos(pkin(9)) + pkin(2);
	t9 = qJ(1) + qJ(2);
	t7 = sin(t9);
	t8 = cos(t9);
	t13 = t16 * t8 + t17 * t7;
	t12 = t16 * t7 - t17 * t8;
	t1 = [0, 0, 0, 0, 0; cos(qJ(1)) * pkin(1) + t13, t13, -t8, 0, 0; sin(qJ(1)) * pkin(1) + t12, t12, -t7, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (71->14), mult. (44->13), div. (0->0), fcn. (46->7), ass. (0->13)
	t22 = r_i_i_C(3) + pkin(7) + qJ(3);
	t12 = pkin(9) + qJ(4);
	t8 = sin(t12);
	t9 = cos(t12);
	t21 = t9 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t20 = cos(pkin(9)) * pkin(3) + pkin(2) + t21;
	t17 = r_i_i_C(1) * t8 + r_i_i_C(2) * t9;
	t13 = qJ(1) + qJ(2);
	t10 = sin(t13);
	t11 = cos(t13);
	t16 = t20 * t10 - t22 * t11;
	t15 = t22 * t10 + t20 * t11;
	t1 = [0, 0, 0, t21, 0; cos(qJ(1)) * pkin(1) + t15, t15, -t11, -t17 * t10, 0; sin(qJ(1)) * pkin(1) + t16, t16, -t10, t17 * t11, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (124->21), mult. (63->22), div. (0->0), fcn. (65->9), ass. (0->19)
	t30 = r_i_i_C(3) + pkin(8) + pkin(7) + qJ(3);
	t19 = pkin(9) + qJ(4);
	t15 = qJ(5) + t19;
	t11 = sin(t15);
	t12 = cos(t15);
	t24 = t12 * r_i_i_C(1) - t11 * r_i_i_C(2);
	t20 = qJ(1) + qJ(2);
	t17 = cos(t20);
	t25 = t12 * t17;
	t26 = t11 * t17;
	t29 = r_i_i_C(1) * t26 + r_i_i_C(2) * t25;
	t28 = pkin(4) * sin(t19);
	t23 = -r_i_i_C(1) * t11 - r_i_i_C(2) * t12;
	t16 = sin(t20);
	t10 = pkin(4) * cos(t19);
	t3 = t10 + cos(pkin(9)) * pkin(3) + pkin(2);
	t22 = -t30 * t17 + (t3 + t24) * t16;
	t21 = r_i_i_C(1) * t25 - r_i_i_C(2) * t26 + t30 * t16 + t17 * t3;
	t1 = [0, 0, 0, t10 + t24, t24; cos(qJ(1)) * pkin(1) + t21, t21, -t17, (t23 - t28) * t16, t23 * t16; sin(qJ(1)) * pkin(1) + t22, t22, -t16, t17 * t28 + t29, t29;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end