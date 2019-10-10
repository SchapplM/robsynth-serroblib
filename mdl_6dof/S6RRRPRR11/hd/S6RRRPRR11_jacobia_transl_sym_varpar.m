% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:10
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR11_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR11_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR11_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
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
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (19->14), mult. (46->24), div. (0->0), fcn. (54->6), ass. (0->16)
	t7 = sin(qJ(2));
	t8 = sin(qJ(1));
	t15 = t8 * t7;
	t9 = cos(qJ(2));
	t14 = t8 * t9;
	t10 = cos(qJ(1));
	t13 = t10 * t7;
	t12 = t10 * t9;
	t5 = sin(pkin(6));
	t11 = (pkin(8) + r_i_i_C(3)) * t5;
	t6 = cos(pkin(6));
	t4 = -t6 * t15 + t12;
	t3 = -t6 * t14 - t13;
	t2 = -t6 * t13 - t14;
	t1 = -t6 * t12 + t15;
	t16 = [-t8 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t10 * t11, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0, 0, 0; t10 * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t11 * t8, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0, 0, 0; 0, (r_i_i_C(1) * t9 - r_i_i_C(2) * t7) * t5, 0, 0, 0, 0;];
	Ja_transl = t16;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (61->29), mult. (150->51), div. (0->0), fcn. (184->8), ass. (0->26)
	t27 = -r_i_i_C(3) - pkin(9);
	t13 = sin(qJ(1));
	t9 = sin(pkin(6));
	t26 = t13 * t9;
	t14 = cos(qJ(3));
	t25 = t14 * t9;
	t16 = cos(qJ(1));
	t24 = t16 * t9;
	t12 = sin(qJ(2));
	t23 = t12 * t13;
	t22 = t12 * t16;
	t15 = cos(qJ(2));
	t21 = t13 * t15;
	t20 = t15 * t16;
	t11 = sin(qJ(3));
	t10 = cos(pkin(6));
	t4 = t10 * t22 + t21;
	t19 = t11 * t24 - t14 * t4;
	t18 = t14 * r_i_i_C(1) - t11 * r_i_i_C(2) + pkin(2);
	t17 = t4 * t11 + t14 * t24;
	t6 = -t10 * t23 + t20;
	t5 = t10 * t21 + t22;
	t3 = -t10 * t20 + t23;
	t2 = t11 * t26 + t14 * t6;
	t1 = -t11 * t6 + t13 * t25;
	t7 = [-t13 * pkin(1) - t4 * pkin(2) + pkin(8) * t24 + t19 * r_i_i_C(1) + t17 * r_i_i_C(2) + t27 * t3, -t18 * t5 - t27 * t6, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; pkin(1) * t16 + t6 * pkin(2) + pkin(8) * t26 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t27 * t5, -t18 * t3 - t27 * t4, -t17 * r_i_i_C(1) + t19 * r_i_i_C(2), 0, 0, 0; 0, (-t27 * t12 + t18 * t15) * t9, (-t11 * t12 * t9 + t10 * t14) * r_i_i_C(1) + (-t10 * t11 - t12 * t25) * r_i_i_C(2), 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (104->31), mult. (255->53), div. (0->0), fcn. (319->8), ass. (0->27)
	t18 = sin(qJ(3));
	t21 = cos(qJ(3));
	t29 = r_i_i_C(3) + qJ(4);
	t34 = pkin(3) + r_i_i_C(1);
	t35 = t29 * t18 + t34 * t21 + pkin(2);
	t33 = pkin(9) + r_i_i_C(2);
	t17 = sin(pkin(6));
	t20 = sin(qJ(1));
	t32 = t17 * t20;
	t31 = t17 * t21;
	t23 = cos(qJ(1));
	t30 = t17 * t23;
	t28 = cos(pkin(6));
	t19 = sin(qJ(2));
	t22 = cos(qJ(2));
	t25 = t23 * t28;
	t10 = t19 * t25 + t20 * t22;
	t27 = t10 * t21 - t18 * t30;
	t26 = t20 * t28;
	t1 = t10 * t18 + t21 * t30;
	t12 = -t19 * t26 + t23 * t22;
	t11 = t23 * t19 + t22 * t26;
	t9 = t20 * t19 - t22 * t25;
	t7 = t17 * t19 * t18 - t28 * t21;
	t6 = t12 * t21 + t18 * t32;
	t5 = t12 * t18 - t20 * t31;
	t2 = [-t20 * pkin(1) - t10 * pkin(2) + pkin(8) * t30 - t29 * t1 - t34 * t27 - t33 * t9, -t11 * t35 + t33 * t12, t29 * t6 - t34 * t5, t5, 0, 0; t23 * pkin(1) + t12 * pkin(2) + pkin(8) * t32 + t33 * t11 + t29 * t5 + t34 * t6, t33 * t10 - t35 * t9, -t34 * t1 + t29 * t27, t1, 0, 0; 0, (t33 * t19 + t35 * t22) * t17, t29 * (t28 * t18 + t19 * t31) - t34 * t7, t7, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (199->45), mult. (502->77), div. (0->0), fcn. (645->10), ass. (0->33)
	t22 = sin(qJ(3));
	t26 = cos(qJ(3));
	t21 = sin(qJ(5));
	t25 = cos(qJ(5));
	t39 = pkin(3) + pkin(4);
	t29 = t25 * r_i_i_C(1) - t21 * r_i_i_C(2) + t39;
	t30 = t21 * r_i_i_C(1) + t25 * r_i_i_C(2) + qJ(4);
	t40 = t30 * t22 + t29 * t26 + pkin(2);
	t38 = cos(qJ(1));
	t20 = sin(pkin(6));
	t24 = sin(qJ(1));
	t37 = t20 * t24;
	t36 = t20 * t26;
	t35 = cos(pkin(6));
	t34 = r_i_i_C(3) + pkin(10) - pkin(9);
	t33 = t20 * t38;
	t23 = sin(qJ(2));
	t27 = cos(qJ(2));
	t31 = t35 * t38;
	t12 = t23 * t31 + t24 * t27;
	t4 = t12 * t26 - t22 * t33;
	t32 = t24 * t35;
	t3 = t12 * t22 + t26 * t33;
	t14 = -t23 * t32 + t38 * t27;
	t13 = t38 * t23 + t27 * t32;
	t11 = t24 * t23 - t27 * t31;
	t10 = t35 * t22 + t23 * t36;
	t9 = t20 * t23 * t22 - t35 * t26;
	t8 = t14 * t26 + t22 * t37;
	t7 = t14 * t22 - t24 * t36;
	t2 = t7 * t21 + t8 * t25;
	t1 = -t8 * t21 + t7 * t25;
	t5 = [-t24 * pkin(1) - t12 * pkin(2) + pkin(8) * t33 + t34 * t11 - t29 * t4 - t30 * t3, -t13 * t40 - t34 * t14, -t29 * t7 + t30 * t8, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t38 * pkin(1) + t14 * pkin(2) + pkin(8) * t37 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(4) - t34 * t13 + t39 * t8, -t11 * t40 - t34 * t12, -t29 * t3 + t30 * t4, t3, (-t4 * t21 + t3 * t25) * r_i_i_C(1) + (-t3 * t21 - t4 * t25) * r_i_i_C(2), 0; 0, (-t34 * t23 + t40 * t27) * t20, t30 * t10 - t29 * t9, t9, (-t10 * t21 + t9 * t25) * r_i_i_C(1) + (-t10 * t25 - t9 * t21) * r_i_i_C(2), 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:19
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (409->67), mult. (1057->107), div. (0->0), fcn. (1379->12), ass. (0->47)
	t43 = sin(qJ(5));
	t44 = sin(qJ(3));
	t48 = cos(qJ(5));
	t49 = cos(qJ(3));
	t42 = sin(qJ(6));
	t47 = cos(qJ(6));
	t55 = t47 * r_i_i_C(1) - t42 * r_i_i_C(2) + pkin(5);
	t72 = r_i_i_C(3) + pkin(11);
	t81 = t72 * (t43 * t49 - t44 * t48) + (t43 * t44 + t48 * t49) * t55;
	t45 = sin(qJ(2));
	t46 = sin(qJ(1));
	t50 = cos(qJ(2));
	t67 = cos(pkin(6));
	t71 = cos(qJ(1));
	t60 = t67 * t71;
	t33 = t45 * t60 + t46 * t50;
	t41 = sin(pkin(6));
	t63 = t41 * t71;
	t22 = t33 * t44 + t49 * t63;
	t23 = t33 * t49 - t44 * t63;
	t80 = t22 * t48 - t23 * t43;
	t6 = t22 * t43 + t23 * t48;
	t74 = pkin(3) + pkin(4);
	t77 = t44 * qJ(4) + t74 * t49 + pkin(2);
	t79 = -t77 - t81;
	t61 = t46 * t67;
	t35 = -t45 * t61 + t71 * t50;
	t69 = t41 * t49;
	t26 = t35 * t44 - t46 * t69;
	t70 = t41 * t46;
	t27 = t35 * t49 + t44 * t70;
	t11 = -t26 * t48 + t27 * t43;
	t12 = t26 * t43 + t27 * t48;
	t78 = -t55 * t11 + t72 * t12;
	t30 = t41 * t45 * t44 - t67 * t49;
	t31 = t67 * t44 + t45 * t69;
	t20 = t30 * t43 + t31 * t48;
	t76 = t55 * (t30 * t48 - t31 * t43) + t72 * t20;
	t75 = t55 * t80 + t72 * t6;
	t73 = pkin(9) - pkin(10);
	t68 = t41 * t50;
	t54 = t42 * r_i_i_C(1) + t47 * r_i_i_C(2) - t73;
	t34 = t71 * t45 + t50 * t61;
	t32 = t46 * t45 - t50 * t60;
	t2 = t12 * t47 - t34 * t42;
	t1 = -t12 * t42 - t34 * t47;
	t3 = [-t46 * pkin(1) - t33 * pkin(2) + pkin(8) * t63 - t22 * qJ(4) - t74 * t23 + t54 * t32 - t55 * t6 + t72 * t80, t79 * t34 - t35 * t54, t27 * qJ(4) - t74 * t26 - t78, t26, t78, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t71 * pkin(1) + t35 * pkin(2) + t12 * pkin(5) + pkin(8) * t70 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t26 * qJ(4) + t72 * t11 + t74 * t27 + t73 * t34, t79 * t32 - t33 * t54, t23 * qJ(4) - t74 * t22 - t75, t22, t75, (-t32 * t47 - t6 * t42) * r_i_i_C(1) + (t32 * t42 - t6 * t47) * r_i_i_C(2); 0, (-t54 * t45 + t77 * t50) * t41 + t81 * t68, t31 * qJ(4) - t74 * t30 - t76, t30, t76, (-t20 * t42 + t47 * t68) * r_i_i_C(1) + (-t20 * t47 - t42 * t68) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end