% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:54
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:54:01
	% EndTime: 2019-10-10 10:54:01
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:54:01
	% EndTime: 2019-10-10 10:54:01
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:54:01
	% EndTime: 2019-10-10 10:54:01
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
	% StartTime: 2019-10-10 10:54:01
	% EndTime: 2019-10-10 10:54:01
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (29->10), mult. (31->12), div. (0->0), fcn. (33->6), ass. (0->10)
	t5 = qJ(2) + pkin(11);
	t2 = sin(t5);
	t3 = cos(t5);
	t14 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2) + cos(qJ(2)) * pkin(2);
	t13 = r_i_i_C(3) + qJ(3) + pkin(7);
	t11 = pkin(1) + t14;
	t10 = -sin(qJ(2)) * pkin(2) - r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t11 * t8 + t13 * t9, t10 * t9, t8, 0, 0, 0; t11 * t9 + t13 * t8, t10 * t8, -t9, 0, 0, 0; 0, t14, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:54:01
	% EndTime: 2019-10-10 10:54:01
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (66->13), mult. (46->16), div. (0->0), fcn. (48->8), ass. (0->13)
	t10 = qJ(2) + pkin(11);
	t7 = qJ(4) + t10;
	t5 = sin(t7);
	t6 = cos(t7);
	t16 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t20 = t16 + pkin(3) * cos(t10) + cos(qJ(2)) * pkin(2);
	t18 = r_i_i_C(3) + pkin(8) + qJ(3) + pkin(7);
	t15 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
	t14 = pkin(1) + t20;
	t13 = t15 - pkin(3) * sin(t10) - sin(qJ(2)) * pkin(2);
	t12 = cos(qJ(1));
	t11 = sin(qJ(1));
	t1 = [-t14 * t11 + t18 * t12, t13 * t12, t11, t15 * t12, 0, 0; t18 * t11 + t14 * t12, t13 * t11, -t12, t15 * t11, 0, 0; 0, t20, 0, t16, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:54:01
	% EndTime: 2019-10-10 10:54:02
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (159->32), mult. (126->41), div. (0->0), fcn. (136->10), ass. (0->30)
	t23 = qJ(2) + pkin(11);
	t20 = qJ(4) + t23;
	t18 = sin(t20);
	t19 = cos(t20);
	t24 = sin(qJ(5));
	t42 = r_i_i_C(2) * t24;
	t48 = pkin(9) + r_i_i_C(3);
	t49 = t18 * t42 + t19 * t48;
	t46 = t19 * pkin(4) + t48 * t18;
	t35 = pkin(3) * cos(t23) + cos(qJ(2)) * pkin(2);
	t45 = pkin(1) + t35 + t46;
	t26 = cos(qJ(5));
	t43 = r_i_i_C(1) * t26;
	t27 = cos(qJ(1));
	t39 = t24 * t27;
	t25 = sin(qJ(1));
	t38 = t25 * t24;
	t37 = t25 * t26;
	t36 = t26 * t27;
	t34 = t49 * t25;
	t32 = t49 * t27;
	t30 = (-pkin(4) - t43) * t18;
	t29 = -pkin(3) * sin(t23) - sin(qJ(2)) * pkin(2) + t30;
	t28 = (-t42 + t43) * t19 + t46;
	t22 = -pkin(8) - qJ(3) - pkin(7);
	t4 = t19 * t36 + t38;
	t3 = -t19 * t39 + t37;
	t2 = -t19 * t37 + t39;
	t1 = t19 * t38 + t36;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t22 * t27 - t45 * t25, t29 * t27 + t32, t25, t27 * t30 + t32, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t25 * t22 + t45 * t27, t29 * t25 + t34, -t27, t25 * t30 + t34, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, t28 + t35, 0, t28, (-r_i_i_C(1) * t24 - r_i_i_C(2) * t26) * t18, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:54:01
	% EndTime: 2019-10-10 10:54:02
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (247->42), mult. (170->54), div. (0->0), fcn. (184->12), ass. (0->40)
	t27 = qJ(2) + pkin(11);
	t22 = qJ(4) + t27;
	t19 = sin(t22);
	t20 = cos(t22);
	t31 = cos(qJ(5));
	t21 = t31 * pkin(5) + pkin(4);
	t33 = -pkin(10) - pkin(9);
	t57 = t20 * t21 + (r_i_i_C(3) - t33) * t19;
	t41 = pkin(3) * cos(t27) + cos(qJ(2)) * pkin(2);
	t56 = pkin(1) + t41 + t57;
	t28 = qJ(5) + qJ(6);
	t24 = cos(t28);
	t32 = cos(qJ(1));
	t43 = t32 * t24;
	t23 = sin(t28);
	t30 = sin(qJ(1));
	t46 = t30 * t23;
	t5 = t20 * t46 + t43;
	t44 = t32 * t23;
	t45 = t30 * t24;
	t6 = -t20 * t45 + t44;
	t55 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t7 = -t20 * t44 + t45;
	t8 = t20 * t43 + t46;
	t54 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t29 = sin(qJ(5));
	t53 = pkin(5) * t29;
	t52 = r_i_i_C(1) * t24;
	t51 = r_i_i_C(2) * t23;
	t40 = t19 * t51;
	t48 = t20 * t30;
	t49 = r_i_i_C(3) * t48 + t30 * t40;
	t47 = t20 * t32;
	t42 = r_i_i_C(3) * t47 + t32 * t40;
	t39 = pkin(8) + qJ(3) + pkin(7) + t53;
	t37 = -r_i_i_C(1) * t23 - r_i_i_C(2) * t24;
	t36 = -t20 * t33 + (-t21 - t52) * t19;
	t35 = (-t51 + t52) * t20 + t57;
	t34 = -pkin(3) * sin(t27) - sin(qJ(2)) * pkin(2) + t36;
	t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t56 * t30 + t39 * t32, t34 * t32 + t42, t30, t36 * t32 + t42, (-t29 * t47 + t30 * t31) * pkin(5) + t54, t54; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t39 * t30 + t56 * t32, t34 * t30 + t49, -t32, t36 * t30 + t49, (-t29 * t48 - t31 * t32) * pkin(5) + t55, t55; 0, t35 + t41, 0, t35, (t37 - t53) * t19, t37 * t19;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end