% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRRR10_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR10_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR10_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:11:37
	% EndTime: 2019-12-31 19:11:37
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:11:37
	% EndTime: 2019-12-31 19:11:37
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:11:37
	% EndTime: 2019-12-31 19:11:37
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (9->6), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(9)) - r_i_i_C(2) * sin(pkin(9)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t3 * t5 + t4 * t6, t3, 0, 0, 0; t3 * t6 + t4 * t5, -t4, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:11:37
	% EndTime: 2019-12-31 19:11:37
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (26->9), mult. (26->11), div. (0->0), fcn. (28->5), ass. (0->10)
	t11 = r_i_i_C(3) + pkin(6) + qJ(2);
	t4 = pkin(9) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t10 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t9 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t8 = cos(pkin(9)) * pkin(2) + pkin(1) + t10;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [t11 * t7 - t8 * t6, t6, t9 * t7, 0, 0; t11 * t6 + t8 * t7, -t7, t9 * t6, 0, 0; 0, 0, t10, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:11:37
	% EndTime: 2019-12-31 19:11:37
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (71->24), mult. (85->33), div. (0->0), fcn. (95->7), ass. (0->22)
	t22 = pkin(7) + r_i_i_C(3);
	t8 = pkin(9) + qJ(3);
	t6 = sin(t8);
	t17 = t22 * t6;
	t7 = cos(t8);
	t23 = t17 + t7 * pkin(3) + cos(pkin(9)) * pkin(2) + pkin(1);
	t10 = sin(qJ(4));
	t13 = cos(qJ(1));
	t21 = t10 * t13;
	t11 = sin(qJ(1));
	t20 = t11 * t10;
	t12 = cos(qJ(4));
	t19 = t11 * t12;
	t18 = t12 * t13;
	t15 = r_i_i_C(1) * t12 - r_i_i_C(2) * t10 + pkin(3);
	t14 = -t15 * t6 + t22 * t7;
	t9 = -pkin(6) - qJ(2);
	t4 = t7 * t18 + t20;
	t3 = -t7 * t21 + t19;
	t2 = -t7 * t19 + t21;
	t1 = t7 * t20 + t18;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t23 * t11 - t13 * t9, t11, t14 * t13, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t11 * t9 + t23 * t13, -t13, t14 * t11, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, 0, t15 * t7 + t17, (-r_i_i_C(1) * t10 - r_i_i_C(2) * t12) * t6, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:11:37
	% EndTime: 2019-12-31 19:11:37
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (138->31), mult. (126->43), div. (0->0), fcn. (140->9), ass. (0->31)
	t20 = cos(qJ(4));
	t10 = pkin(4) * t20 + pkin(3);
	t15 = pkin(9) + qJ(3);
	t12 = cos(t15);
	t11 = sin(t15);
	t34 = r_i_i_C(3) + pkin(8) + pkin(7);
	t28 = t34 * t11;
	t38 = t28 + t12 * t10 + cos(pkin(9)) * pkin(2) + pkin(1);
	t16 = qJ(4) + qJ(5);
	t14 = cos(t16);
	t21 = cos(qJ(1));
	t29 = t14 * t21;
	t13 = sin(t16);
	t19 = sin(qJ(1));
	t32 = t13 * t19;
	t5 = t12 * t32 + t29;
	t30 = t14 * t19;
	t31 = t13 * t21;
	t6 = -t12 * t30 + t31;
	t37 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t7 = -t12 * t31 + t30;
	t8 = t12 * t29 + t32;
	t36 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t18 = sin(qJ(4));
	t35 = pkin(4) * t18;
	t33 = t12 * t18;
	t27 = pkin(6) + qJ(2) + t35;
	t25 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
	t24 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t10;
	t23 = -t24 * t11 + t34 * t12;
	t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t38 * t19 + t27 * t21, t19, t23 * t21, (t19 * t20 - t21 * t33) * pkin(4) + t36, t36; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t27 * t19 + t38 * t21, -t21, t23 * t19, (-t19 * t33 - t20 * t21) * pkin(4) + t37, t37; 0, 0, t24 * t12 + t28, (t25 - t35) * t11, t25 * t11;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end