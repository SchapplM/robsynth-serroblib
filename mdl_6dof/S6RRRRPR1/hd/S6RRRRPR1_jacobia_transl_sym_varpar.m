% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:35
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:06
	% EndTime: 2019-10-10 12:35:06
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:06
	% EndTime: 2019-10-10 12:35:06
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:06
	% EndTime: 2019-10-10 12:35:06
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (11->6), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->9)
	t8 = pkin(7) + r_i_i_C(3);
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t2 * t5 + t4 * t8, t6 * t4, 0, 0, 0, 0; t2 * t8 + t4 * t5, t6 * t2, 0, 0, 0, 0; 0, t7, 0, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:06
	% EndTime: 2019-10-10 12:35:06
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (37->9), mult. (41->14), div. (0->0), fcn. (41->6), ass. (0->12)
	t6 = qJ(2) + qJ(3);
	t3 = sin(t6);
	t4 = cos(t6);
	t14 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t17 = t14 + cos(qJ(2)) * pkin(2);
	t15 = r_i_i_C(3) + pkin(8) + pkin(7);
	t13 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t12 = pkin(1) + t17;
	t11 = -sin(qJ(2)) * pkin(2) + t13;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t12 * t8 + t15 * t9, t11 * t9, t13 * t9, 0, 0, 0; t12 * t9 + t15 * t8, t11 * t8, t13 * t8, 0, 0, 0; 0, t17, t14, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:06
	% EndTime: 2019-10-10 12:35:06
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (86->12), mult. (61->18), div. (0->0), fcn. (61->8), ass. (0->15)
	t11 = qJ(2) + qJ(3);
	t8 = qJ(4) + t11;
	t5 = sin(t8);
	t6 = cos(t8);
	t19 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t18 = t19 + pkin(3) * cos(t11);
	t23 = t18 + cos(qJ(2)) * pkin(2);
	t17 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
	t14 = t17 - pkin(3) * sin(t11);
	t20 = r_i_i_C(3) + pkin(9) + pkin(8) + pkin(7);
	t16 = pkin(1) + t23;
	t15 = -sin(qJ(2)) * pkin(2) + t14;
	t13 = cos(qJ(1));
	t12 = sin(qJ(1));
	t1 = [-t16 * t12 + t20 * t13, t15 * t13, t14 * t13, t17 * t13, 0, 0; t20 * t12 + t16 * t13, t15 * t12, t14 * t12, t17 * t12, 0, 0; 0, t23, t18, t19, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:06
	% EndTime: 2019-10-10 12:35:06
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (144->16), mult. (74->20), div. (0->0), fcn. (76->10), ass. (0->16)
	t15 = qJ(2) + qJ(3);
	t13 = qJ(4) + t15;
	t9 = pkin(11) + t13;
	t5 = sin(t9);
	t6 = cos(t9);
	t24 = t6 * r_i_i_C(1) - r_i_i_C(2) * t5 + pkin(4) * cos(t13);
	t22 = t24 + pkin(3) * cos(t15);
	t28 = cos(qJ(2)) * pkin(2) + t22;
	t18 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6 - pkin(4) * sin(t13);
	t19 = -pkin(3) * sin(t15) + t18;
	t25 = r_i_i_C(3) + qJ(5) + pkin(9) + pkin(8) + pkin(7);
	t21 = pkin(1) + t28;
	t20 = -sin(qJ(2)) * pkin(2) + t19;
	t17 = cos(qJ(1));
	t16 = sin(qJ(1));
	t1 = [-t21 * t16 + t25 * t17, t20 * t17, t19 * t17, t18 * t17, t16, 0; t25 * t16 + t21 * t17, t20 * t16, t19 * t16, t18 * t16, -t17, 0; 0, t28, t22, t24, 0, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:06
	% EndTime: 2019-10-10 12:35:07
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (303->37), mult. (175->45), div. (0->0), fcn. (185->12), ass. (0->34)
	t28 = qJ(2) + qJ(3);
	t26 = qJ(4) + t28;
	t22 = pkin(11) + t26;
	t18 = sin(t22);
	t19 = cos(t22);
	t29 = sin(qJ(6));
	t49 = r_i_i_C(2) * t29;
	t57 = pkin(10) + r_i_i_C(3);
	t58 = t18 * t49 + t19 * t57;
	t31 = cos(qJ(6));
	t50 = r_i_i_C(1) * t31;
	t35 = (-pkin(5) - t50) * t18 - pkin(4) * sin(t26);
	t55 = t57 * t18 + t19 * pkin(5) + pkin(4) * cos(t26);
	t36 = -pkin(3) * sin(t28) + t35;
	t21 = pkin(3) * cos(t28);
	t27 = cos(qJ(2)) * pkin(2);
	t53 = pkin(1) + t21 + t27 + t55;
	t32 = cos(qJ(1));
	t46 = t29 * t32;
	t30 = sin(qJ(1));
	t45 = t30 * t29;
	t44 = t30 * t31;
	t43 = t31 * t32;
	t41 = t58 * t30;
	t40 = t58 * t32;
	t37 = -sin(qJ(2)) * pkin(2) + t36;
	t34 = (-t49 + t50) * t19 + t55;
	t33 = t21 + t34;
	t25 = -qJ(5) - pkin(9) - pkin(8) - pkin(7);
	t4 = t19 * t43 + t45;
	t3 = -t19 * t46 + t44;
	t2 = -t19 * t44 + t46;
	t1 = t19 * t45 + t43;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t25 * t32 - t53 * t30, t37 * t32 + t40, t36 * t32 + t40, t35 * t32 + t40, t30, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t30 * t25 + t53 * t32, t37 * t30 + t41, t36 * t30 + t41, t35 * t30 + t41, -t32, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, t27 + t33, t33, t34, 0, (-r_i_i_C(1) * t29 - r_i_i_C(2) * t31) * t18;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end