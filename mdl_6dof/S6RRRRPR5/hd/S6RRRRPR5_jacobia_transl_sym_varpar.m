% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:38
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:38:37
	% EndTime: 2019-10-10 12:38:37
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:38:37
	% EndTime: 2019-10-10 12:38:37
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
	% StartTime: 2019-10-10 12:38:37
	% EndTime: 2019-10-10 12:38:37
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
	% StartTime: 2019-10-10 12:38:37
	% EndTime: 2019-10-10 12:38:37
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
	% StartTime: 2019-10-10 12:38:37
	% EndTime: 2019-10-10 12:38:37
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (100->28), mult. (121->39), div. (0->0), fcn. (129->8), ass. (0->29)
	t19 = qJ(2) + qJ(3);
	t16 = sin(t19);
	t17 = cos(t19);
	t20 = sin(qJ(4));
	t39 = r_i_i_C(2) * t20;
	t45 = pkin(9) + r_i_i_C(3);
	t46 = t16 * t39 + t17 * t45;
	t43 = t17 * pkin(3) + t45 * t16;
	t18 = cos(qJ(2)) * pkin(2);
	t42 = pkin(1) + t18 + t43;
	t23 = cos(qJ(4));
	t40 = r_i_i_C(1) * t23;
	t22 = sin(qJ(1));
	t36 = t20 * t22;
	t24 = cos(qJ(1));
	t35 = t20 * t24;
	t34 = t22 * t23;
	t33 = t23 * t24;
	t32 = t46 * t22;
	t30 = t46 * t24;
	t28 = (-pkin(3) - t40) * t16;
	t27 = (-t39 + t40) * t17 + t43;
	t26 = -sin(qJ(2)) * pkin(2) + t28;
	t25 = -pkin(8) - pkin(7);
	t4 = t17 * t33 + t36;
	t3 = -t17 * t35 + t34;
	t2 = -t17 * t34 + t35;
	t1 = t17 * t36 + t33;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t42 * t22 - t24 * t25, t26 * t24 + t30, t24 * t28 + t30, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t22 * t25 + t42 * t24, t26 * t22 + t32, t22 * t28 + t32, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, t18 + t27, t27, (-r_i_i_C(1) * t20 - r_i_i_C(2) * t23) * t16, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:38:37
	% EndTime: 2019-10-10 12:38:37
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (157->29), mult. (196->39), div. (0->0), fcn. (217->8), ass. (0->30)
	t23 = sin(qJ(4));
	t26 = cos(qJ(4));
	t34 = r_i_i_C(3) + qJ(5);
	t45 = pkin(4) + r_i_i_C(1);
	t51 = t34 * t23 + t45 * t26;
	t50 = pkin(9) + r_i_i_C(2);
	t22 = qJ(2) + qJ(3);
	t20 = cos(t22);
	t48 = t20 * t50;
	t19 = sin(t22);
	t47 = t20 * pkin(3) + t50 * t19;
	t21 = cos(qJ(2)) * pkin(2);
	t46 = pkin(1) + t21 + t47;
	t25 = sin(qJ(1));
	t43 = t25 * t48;
	t38 = t25 * t23;
	t37 = t25 * t26;
	t27 = cos(qJ(1));
	t36 = t27 * t23;
	t35 = t27 * t26;
	t33 = t27 * t48;
	t31 = t51 * t20 + t47;
	t30 = (-pkin(3) - t51) * t19;
	t29 = -sin(qJ(2)) * pkin(2) + t30;
	t28 = -pkin(8) - pkin(7);
	t4 = t20 * t35 + t38;
	t3 = t20 * t36 - t37;
	t2 = t20 * t37 - t36;
	t1 = t20 * t38 + t35;
	t5 = [-t34 * t1 - t45 * t2 - t46 * t25 - t27 * t28, t29 * t27 + t33, t27 * t30 + t33, -t45 * t3 + t34 * t4, t3, 0; -t25 * t28 + t46 * t27 + t34 * t3 + t45 * t4, t29 * t25 + t43, t25 * t30 + t43, -t45 * t1 + t34 * t2, t1, 0; 0, t21 + t31, t31, (-t45 * t23 + t34 * t26) * t19, t19 * t23, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:38:37
	% EndTime: 2019-10-10 12:38:37
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (280->52), mult. (380->72), div. (0->0), fcn. (449->10), ass. (0->43)
	t31 = sin(qJ(4));
	t35 = cos(qJ(4));
	t62 = pkin(4) + pkin(5);
	t67 = -qJ(5) * t31 - t62 * t35;
	t29 = qJ(2) + qJ(3);
	t27 = cos(t29);
	t26 = sin(t29);
	t30 = sin(qJ(6));
	t34 = cos(qJ(6));
	t45 = t30 * t31 + t34 * t35;
	t64 = t45 * t26;
	t46 = t30 * t35 - t31 * t34;
	t9 = t46 * t26;
	t66 = pkin(9) * t27 - r_i_i_C(1) * t64 + r_i_i_C(2) * t9;
	t61 = -pkin(10) - r_i_i_C(3);
	t65 = t27 * pkin(3) + (pkin(9) + t61) * t26;
	t28 = cos(qJ(2)) * pkin(2);
	t63 = pkin(1) + t28 + t65;
	t36 = cos(qJ(1));
	t55 = t36 * t35;
	t33 = sin(qJ(1));
	t58 = t33 * t31;
	t12 = t27 * t58 + t55;
	t11 = t12 * t34;
	t57 = t33 * t35;
	t56 = t36 * t31;
	t52 = t66 * t33;
	t51 = t66 * t36;
	t49 = -t30 * r_i_i_C(1) - qJ(5);
	t48 = t30 * r_i_i_C(2) - t62;
	t47 = -t9 * r_i_i_C(1) - r_i_i_C(2) * t64;
	t44 = t34 * r_i_i_C(2) - t49;
	t43 = -t34 * r_i_i_C(1) + t48;
	t40 = t65 + (r_i_i_C(1) * t45 - r_i_i_C(2) * t46 - t67) * t27;
	t39 = t61 * t27 + (-pkin(3) + t67) * t26;
	t38 = -sin(qJ(2)) * pkin(2) + t39;
	t37 = -pkin(8) - pkin(7);
	t15 = t27 * t55 + t58;
	t14 = t27 * t56 - t57;
	t13 = t27 * t57 - t56;
	t2 = t14 * t30 + t15 * t34;
	t1 = t14 * t34 - t15 * t30;
	t3 = [-t11 * r_i_i_C(2) + t49 * t12 + t43 * t13 - t63 * t33 - t36 * t37, t36 * t38 + t51, t36 * t39 + t51, t43 * t14 + t44 * t15, t14, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t14 * qJ(5) + t62 * t15 - t33 * t37 + t63 * t36, t33 * t38 + t52, t33 * t39 + t52, -t11 * r_i_i_C(1) + t48 * t12 + t44 * t13, t12, (-t13 * t30 + t11) * r_i_i_C(1) + (-t12 * t30 - t13 * t34) * r_i_i_C(2); 0, t28 + t40, t40, (qJ(5) * t35 - t62 * t31) * t26 - t47, t26 * t31, t47;];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end