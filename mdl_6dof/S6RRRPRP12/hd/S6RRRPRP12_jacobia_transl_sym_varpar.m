% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:53
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP12_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP12_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP12_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:53:24
	% EndTime: 2019-10-10 11:53:24
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:53:24
	% EndTime: 2019-10-10 11:53:24
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
	% StartTime: 2019-10-10 11:53:24
	% EndTime: 2019-10-10 11:53:24
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
	% StartTime: 2019-10-10 11:53:24
	% EndTime: 2019-10-10 11:53:24
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
	% StartTime: 2019-10-10 11:53:24
	% EndTime: 2019-10-10 11:53:24
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (104->31), mult. (255->53), div. (0->0), fcn. (319->8), ass. (0->27)
	t16 = sin(qJ(3));
	t19 = cos(qJ(3));
	t27 = r_i_i_C(3) + qJ(4);
	t32 = pkin(3) - r_i_i_C(2);
	t33 = t27 * t16 + t32 * t19 + pkin(2);
	t31 = pkin(9) + r_i_i_C(1);
	t15 = sin(pkin(6));
	t18 = sin(qJ(1));
	t30 = t15 * t18;
	t29 = t15 * t19;
	t21 = cos(qJ(1));
	t28 = t15 * t21;
	t26 = cos(pkin(6));
	t25 = t18 * t26;
	t24 = t21 * t26;
	t17 = sin(qJ(2));
	t20 = cos(qJ(2));
	t10 = t17 * t24 + t18 * t20;
	t1 = t10 * t16 + t19 * t28;
	t23 = -t10 * t19 + t16 * t28;
	t12 = -t17 * t25 + t21 * t20;
	t11 = t21 * t17 + t20 * t25;
	t9 = t18 * t17 - t20 * t24;
	t7 = t15 * t17 * t16 - t26 * t19;
	t6 = t12 * t19 + t16 * t30;
	t5 = t12 * t16 - t18 * t29;
	t2 = [-t18 * pkin(1) - t10 * pkin(2) + pkin(8) * t28 - t27 * t1 + t32 * t23 - t31 * t9, -t11 * t33 + t31 * t12, t27 * t6 - t32 * t5, t5, 0, 0; t21 * pkin(1) + t12 * pkin(2) + pkin(8) * t30 + t31 * t11 + t27 * t5 + t32 * t6, t31 * t10 - t33 * t9, -t32 * t1 - t27 * t23, t1, 0, 0; 0, (t31 * t17 + t33 * t20) * t15, t27 * (t26 * t16 + t17 * t29) - t32 * t7, t7, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:53:24
	% EndTime: 2019-10-10 11:53:24
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (183->45), mult. (459->78), div. (0->0), fcn. (583->10), ass. (0->33)
	t20 = sin(qJ(3));
	t24 = cos(qJ(3));
	t19 = sin(qJ(5));
	t23 = cos(qJ(5));
	t29 = t19 * r_i_i_C(1) + t23 * r_i_i_C(2) + qJ(4);
	t33 = pkin(3) + pkin(10) + r_i_i_C(3);
	t40 = t29 * t20 + t33 * t24 + pkin(2);
	t39 = pkin(4) + pkin(9);
	t38 = cos(qJ(1));
	t18 = sin(pkin(6));
	t22 = sin(qJ(1));
	t37 = t18 * t22;
	t36 = t18 * t24;
	t25 = cos(qJ(2));
	t35 = t18 * t25;
	t34 = cos(pkin(6));
	t32 = t18 * t38;
	t31 = t22 * t34;
	t30 = t34 * t38;
	t28 = t23 * r_i_i_C(1) - t19 * r_i_i_C(2) + t39;
	t21 = sin(qJ(2));
	t12 = t21 * t30 + t22 * t25;
	t3 = t12 * t20 + t24 * t32;
	t27 = t12 * t24 - t20 * t32;
	t14 = -t21 * t31 + t38 * t25;
	t13 = t38 * t21 + t25 * t31;
	t11 = t22 * t21 - t25 * t30;
	t9 = t18 * t21 * t20 - t34 * t24;
	t8 = t14 * t24 + t20 * t37;
	t7 = t14 * t20 - t22 * t36;
	t2 = t13 * t23 + t7 * t19;
	t1 = -t13 * t19 + t7 * t23;
	t4 = [-t22 * pkin(1) - t12 * pkin(2) + pkin(8) * t32 - t28 * t11 - t33 * t27 - t29 * t3, -t13 * t40 + t28 * t14, t29 * t8 - t33 * t7, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t38 * pkin(1) + t14 * pkin(2) + pkin(8) * t37 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(4) + t39 * t13 + t33 * t8, -t11 * t40 + t28 * t12, t29 * t27 - t33 * t3, t3, (-t11 * t19 + t3 * t23) * r_i_i_C(1) + (-t11 * t23 - t3 * t19) * r_i_i_C(2), 0; 0, (t28 * t21 + t40 * t25) * t18, -t33 * t9 + t29 * (t34 * t20 + t21 * t36), t9, (t19 * t35 + t9 * t23) * r_i_i_C(1) + (-t9 * t19 + t23 * t35) * r_i_i_C(2), 0;];
	Ja_transl = t4;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:53:24
	% EndTime: 2019-10-10 11:53:25
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (275->61), mult. (691->99), div. (0->0), fcn. (888->10), ass. (0->40)
	t35 = sin(qJ(2));
	t36 = sin(qJ(1));
	t39 = cos(qJ(2));
	t48 = cos(pkin(6));
	t58 = cos(qJ(1));
	t43 = t48 * t58;
	t26 = t35 * t43 + t36 * t39;
	t34 = sin(qJ(3));
	t38 = cos(qJ(3));
	t32 = sin(pkin(6));
	t46 = t32 * t58;
	t15 = t26 * t34 + t38 * t46;
	t25 = t36 * t35 - t39 * t43;
	t33 = sin(qJ(5));
	t37 = cos(qJ(5));
	t63 = t15 * t33 + t25 * t37;
	t62 = t15 * t37 - t25 * t33;
	t47 = pkin(3) + pkin(10) + r_i_i_C(2);
	t61 = qJ(4) * t34 + t47 * t38 + pkin(2);
	t60 = pkin(4) + pkin(9);
	t59 = pkin(5) + r_i_i_C(1);
	t55 = t32 * t36;
	t54 = t32 * t38;
	t53 = t33 * t34;
	t52 = t33 * t39;
	t51 = t34 * t37;
	t50 = t37 * t39;
	t49 = r_i_i_C(3) + qJ(6);
	t44 = t36 * t48;
	t42 = t26 * t38 - t34 * t46;
	t40 = t59 * t33 - t49 * t37 + qJ(4);
	t28 = -t35 * t44 + t58 * t39;
	t27 = t58 * t35 + t39 * t44;
	t23 = t32 * t35 * t34 - t48 * t38;
	t20 = t28 * t38 + t34 * t55;
	t19 = t28 * t34 - t36 * t54;
	t13 = t23 * t37 + t32 * t52;
	t6 = t19 * t33 + t27 * t37;
	t5 = -t19 * t37 + t27 * t33;
	t1 = [-t36 * pkin(1) - t26 * pkin(2) + pkin(8) * t46 - t15 * qJ(4) - t60 * t25 - t47 * t42 + t49 * t62 - t59 * t63, t49 * (t27 * t51 + t28 * t33) + t60 * t28 + t59 * (-t27 * t53 + t28 * t37) - t61 * t27, -t47 * t19 + t40 * t20, t19, t49 * t6 - t59 * t5, t5; t58 * pkin(1) + t28 * pkin(2) + pkin(8) * t55 + t19 * qJ(4) + t47 * t20 + t60 * t27 + t49 * t5 + t59 * t6, t59 * (-t25 * t53 + t26 * t37) + t49 * (t25 * t51 + t26 * t33) + t60 * t26 - t61 * t25, -t47 * t15 + t40 * t42, t15, t49 * t63 + t59 * t62, -t62; 0, (t59 * (t34 * t52 + t35 * t37) + t49 * (t33 * t35 - t34 * t50) + t60 * t35 + t61 * t39) * t32, -t47 * t23 + t40 * (t48 * t34 + t35 * t54), t23, -t49 * (-t23 * t33 + t32 * t50) + t59 * t13, -t13;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end