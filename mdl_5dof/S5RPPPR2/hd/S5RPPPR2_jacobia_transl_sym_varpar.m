% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPPPR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPPR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:57
	% EndTime: 2020-01-03 11:23:57
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:57
	% EndTime: 2020-01-03 11:23:57
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
	% StartTime: 2020-01-03 11:23:57
	% EndTime: 2020-01-03 11:23:57
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (10->7), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(7)) - r_i_i_C(2) * sin(pkin(7)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t3 * t6 + t4 * t5, -t4, 0, 0, 0; t3 * t5 - t4 * t6, -t3, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:57
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (20->13), mult. (34->12), div. (0->0), fcn. (43->6), ass. (0->9)
	t1 = sin(pkin(8));
	t3 = cos(pkin(8));
	t8 = t1 * r_i_i_C(1) + t3 * r_i_i_C(2) + qJ(2);
	t2 = sin(pkin(7));
	t4 = cos(pkin(7));
	t7 = pkin(1) + (r_i_i_C(3) + qJ(3)) * t2 + (r_i_i_C(1) * t3 - r_i_i_C(2) * t1 + pkin(2)) * t4;
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t9 = [0, 0, -t4, 0, 0; t8 * t5 + t7 * t6, -t6, t5 * t2, 0, 0; t7 * t5 - t8 * t6, -t5, -t6 * t2, 0, 0;];
	Ja_transl = t9;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:57
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (36->22), mult. (73->26), div. (0->0), fcn. (97->8), ass. (0->17)
	t12 = sin(qJ(1));
	t7 = sin(pkin(8));
	t19 = t12 * t7;
	t11 = cos(pkin(7));
	t13 = cos(qJ(1));
	t18 = t11 * t13;
	t10 = cos(pkin(8));
	t17 = t12 * t10;
	t16 = r_i_i_C(3) + qJ(4);
	t6 = sin(pkin(9));
	t9 = cos(pkin(9));
	t15 = r_i_i_C(1) * t9 - r_i_i_C(2) * t6 + pkin(3);
	t8 = sin(pkin(7));
	t14 = t11 * pkin(2) + pkin(1) + (r_i_i_C(1) * t6 + r_i_i_C(2) * t9 + qJ(3)) * t8;
	t3 = t7 * t18 - t17;
	t1 = t10 * t13 + t11 * t19;
	t2 = [0, 0, -t11, t8 * t7, 0; t12 * qJ(2) + t16 * t3 + t15 * (t10 * t18 + t19) + t14 * t13, -t13, t12 * t8, t1, 0; -t13 * qJ(2) + t16 * t1 + t15 * (t11 * t17 - t13 * t7) + t14 * t12, -t12, -t13 * t8, -t3, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:57
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (81->40), mult. (189->59), div. (0->0), fcn. (251->10), ass. (0->30)
	t15 = sin(pkin(8));
	t18 = cos(pkin(8));
	t21 = sin(qJ(1));
	t26 = t21 * t18;
	t19 = cos(pkin(7));
	t23 = cos(qJ(1));
	t29 = t19 * t23;
	t11 = t15 * t29 - t26;
	t20 = sin(qJ(5));
	t22 = cos(qJ(5));
	t28 = t21 * t15;
	t12 = t18 * t29 + t28;
	t14 = sin(pkin(9));
	t17 = cos(pkin(9));
	t16 = sin(pkin(7));
	t25 = t23 * t16;
	t6 = t12 * t17 + t14 * t25;
	t35 = -t11 * t22 + t20 * t6;
	t34 = t11 * t20 + t22 * t6;
	t33 = r_i_i_C(3) + pkin(6);
	t30 = t16 * t15;
	t27 = t21 * t16;
	t24 = t19 * pkin(2) + t16 * qJ(3) + pkin(1);
	t10 = -t15 * t23 + t19 * t26;
	t9 = t18 * t23 + t19 * t28;
	t8 = t16 * t17 * t18 - t14 * t19;
	t4 = t10 * t17 + t14 * t27;
	t2 = t20 * t9 + t22 * t4;
	t1 = -t20 * t4 + t22 * t9;
	t3 = [0, 0, -t19, t30, (-t20 * t8 + t22 * t30) * r_i_i_C(1) + (-t20 * t30 - t22 * t8) * r_i_i_C(2); t34 * r_i_i_C(1) - t35 * r_i_i_C(2) + t6 * pkin(4) + t12 * pkin(3) + t11 * qJ(4) + t21 * qJ(2) + t33 * (t12 * t14 - t17 * t25) + t24 * t23, -t23, t27, t9, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t10 * pkin(3) + t4 * pkin(4) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - qJ(2) * t23 + t9 * qJ(4) + t33 * (t10 * t14 - t17 * t27) + t24 * t21, -t21, -t25, -t11, t35 * r_i_i_C(1) + t34 * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,5);
end