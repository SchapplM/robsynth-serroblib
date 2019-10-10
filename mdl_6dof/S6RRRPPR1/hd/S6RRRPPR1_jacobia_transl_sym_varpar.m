% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:17
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
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
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
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
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
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
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (72->13), mult. (51->16), div. (0->0), fcn. (53->8), ass. (0->13)
	t11 = qJ(2) + qJ(3);
	t7 = pkin(10) + t11;
	t4 = sin(t7);
	t5 = cos(t7);
	t18 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2) + pkin(3) * cos(t11);
	t22 = t18 + cos(qJ(2)) * pkin(2);
	t14 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5 - pkin(3) * sin(t11);
	t19 = r_i_i_C(3) + qJ(4) + pkin(8) + pkin(7);
	t16 = pkin(1) + t22;
	t15 = -sin(qJ(2)) * pkin(2) + t14;
	t13 = cos(qJ(1));
	t12 = sin(qJ(1));
	t1 = [-t16 * t12 + t13 * t19, t15 * t13, t14 * t13, t12, 0, 0; t12 * t19 + t16 * t13, t15 * t12, t14 * t12, -t13, 0, 0; 0, t22, t18, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (153->24), mult. (113->25), div. (0->0), fcn. (122->10), ass. (0->21)
	t45 = r_i_i_C(3) + qJ(5);
	t20 = qJ(2) + qJ(3);
	t16 = pkin(10) + t20;
	t13 = sin(t16);
	t14 = cos(t16);
	t22 = cos(pkin(11));
	t31 = -r_i_i_C(1) * t22 - pkin(4);
	t21 = sin(pkin(11));
	t38 = r_i_i_C(2) * t21;
	t25 = t14 * (-t31 - t38) + pkin(3) * cos(t20) + t45 * t13;
	t44 = cos(qJ(2)) * pkin(2) + t25;
	t43 = t13 * t38 + t45 * t14;
	t26 = t31 * t13 - pkin(3) * sin(t20);
	t41 = pkin(1) + t44;
	t23 = sin(qJ(1));
	t34 = t43 * t23;
	t24 = cos(qJ(1));
	t33 = t43 * t24;
	t28 = t21 * r_i_i_C(1) + t22 * r_i_i_C(2) + pkin(7) + pkin(8) + qJ(4);
	t27 = -sin(qJ(2)) * pkin(2) + t26;
	t1 = [-t41 * t23 + t28 * t24, t24 * t27 + t33, t24 * t26 + t33, t23, t24 * t13, 0; t28 * t23 + t41 * t24, t23 * t27 + t34, t23 * t26 + t34, -t24, t23 * t13, 0; 0, t44, t25, 0, -t14, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (220->39), mult. (145->46), div. (0->0), fcn. (158->12), ass. (0->32)
	t25 = qJ(2) + qJ(3);
	t20 = pkin(10) + t25;
	t14 = sin(t20);
	t15 = cos(t20);
	t24 = pkin(11) + qJ(6);
	t18 = sin(t24);
	t45 = r_i_i_C(2) * t18;
	t51 = r_i_i_C(3) * t15 + t14 * t45;
	t16 = cos(pkin(11)) * pkin(5) + pkin(4);
	t27 = -pkin(9) - qJ(5);
	t50 = (r_i_i_C(3) - t27) * t14 + t15 * t16 + pkin(3) * cos(t25);
	t19 = cos(t24);
	t46 = r_i_i_C(1) * t19;
	t30 = (-t16 - t46) * t14 - t15 * t27 - pkin(3) * sin(t25);
	t22 = cos(qJ(2)) * pkin(2);
	t48 = pkin(1) + t22 + t50;
	t28 = sin(qJ(1));
	t42 = t51 * t28;
	t29 = cos(qJ(1));
	t41 = t51 * t29;
	t40 = t18 * t29;
	t39 = t19 * t29;
	t38 = t28 * t18;
	t37 = t28 * t19;
	t35 = pkin(5) * sin(pkin(11)) + pkin(7) + pkin(8) + qJ(4);
	t32 = -sin(qJ(2)) * pkin(2) + t30;
	t31 = (-t45 + t46) * t15 + t50;
	t4 = t15 * t39 + t38;
	t3 = -t15 * t40 + t37;
	t2 = -t15 * t37 + t40;
	t1 = t15 * t38 + t39;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t48 * t28 + t35 * t29, t29 * t32 + t41, t29 * t30 + t41, t28, t29 * t14, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t35 * t28 + t48 * t29, t28 * t32 + t42, t28 * t30 + t42, -t29, t28 * t14, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, t22 + t31, t31, 0, -t15, (-r_i_i_C(1) * t18 - r_i_i_C(2) * t19) * t14;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end