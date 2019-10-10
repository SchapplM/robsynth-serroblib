% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:53
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:47
	% EndTime: 2019-10-09 21:53:47
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:47
	% EndTime: 2019-10-09 21:53:47
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:47
	% EndTime: 2019-10-09 21:53:47
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(6));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(11));
	t1 = sin(pkin(11));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:47
	% EndTime: 2019-10-09 21:53:47
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (23->15), mult. (60->31), div. (0->0), fcn. (77->8), ass. (0->14)
	t10 = cos(pkin(6));
	t12 = cos(qJ(2));
	t13 = t10 * t12;
	t11 = sin(qJ(2));
	t5 = sin(pkin(12));
	t8 = cos(pkin(12));
	t4 = -t11 * t8 - t12 * t5;
	t3 = t11 * t5 - t12 * t8;
	t9 = cos(pkin(11));
	t7 = sin(pkin(6));
	t6 = sin(pkin(11));
	t2 = t4 * t10;
	t1 = t3 * t10;
	t14 = [0, (t6 * t1 + t9 * t4) * r_i_i_C(1) + (-t6 * t2 + t9 * t3) * r_i_i_C(2) + (-t9 * t11 - t6 * t13) * pkin(2), t6 * t7, 0, 0, 0; 0, (-t9 * t1 + t6 * t4) * r_i_i_C(1) + (t9 * t2 + t6 * t3) * r_i_i_C(2) + (-t6 * t11 + t9 * t13) * pkin(2), -t9 * t7, 0, 0, 0; 1, (t12 * pkin(2) - t3 * r_i_i_C(1) + t4 * r_i_i_C(2)) * t7, t10, 0, 0, 0;];
	Ja_transl = t14;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:47
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (76->27), mult. (197->52), div. (0->0), fcn. (255->10), ass. (0->23)
	t30 = pkin(8) + r_i_i_C(3);
	t15 = sin(pkin(11));
	t16 = sin(pkin(6));
	t29 = t15 * t16;
	t18 = cos(pkin(11));
	t28 = t18 * t16;
	t19 = cos(pkin(6));
	t23 = cos(qJ(2));
	t27 = t19 * t23;
	t14 = sin(pkin(12));
	t17 = cos(pkin(12));
	t21 = sin(qJ(2));
	t25 = t14 * t23 + t21 * t17;
	t10 = t25 * t19;
	t11 = t14 * t21 - t23 * t17;
	t3 = t10 * t18 - t11 * t15;
	t26 = t10 * t15 + t11 * t18;
	t20 = sin(qJ(4));
	t22 = cos(qJ(4));
	t24 = r_i_i_C(1) * t22 - r_i_i_C(2) * t20 + pkin(3);
	t9 = t11 * t19;
	t8 = t25 * t16;
	t1 = [0, -t30 * t26 + (-t15 * t27 - t18 * t21) * pkin(2) + t24 * (t15 * t9 - t18 * t25), t29, (t20 * t26 + t22 * t29) * r_i_i_C(1) + (-t20 * t29 + t22 * t26) * r_i_i_C(2), 0, 0; 0, t30 * t3 + (-t15 * t21 + t18 * t27) * pkin(2) + t24 * (-t15 * t25 - t18 * t9), -t28, (-t20 * t3 - t22 * t28) * r_i_i_C(1) + (t20 * t28 - t22 * t3) * r_i_i_C(2), 0, 0; 1, t30 * t8 + (pkin(2) * t23 - t11 * t24) * t16, t19, (t19 * t22 - t20 * t8) * r_i_i_C(1) + (-t19 * t20 - t22 * t8) * r_i_i_C(2), 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:47
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (148->35), mult. (296->62), div. (0->0), fcn. (381->12), ass. (0->29)
	t23 = qJ(4) + qJ(5);
	t21 = sin(t23);
	t22 = cos(t23);
	t26 = sin(pkin(6));
	t28 = cos(pkin(11));
	t40 = t28 * t26;
	t29 = cos(pkin(6));
	t24 = sin(pkin(12));
	t27 = cos(pkin(12));
	t31 = sin(qJ(2));
	t33 = cos(qJ(2));
	t36 = t33 * t24 + t31 * t27;
	t16 = t36 * t29;
	t17 = t31 * t24 - t33 * t27;
	t25 = sin(pkin(11));
	t7 = t28 * t16 - t25 * t17;
	t44 = (-t7 * t21 - t22 * t40) * r_i_i_C(1) + (t21 * t40 - t7 * t22) * r_i_i_C(2);
	t37 = t25 * t16 + t28 * t17;
	t41 = t25 * t26;
	t43 = (t21 * t37 + t22 * t41) * r_i_i_C(1) + (-t21 * t41 + t22 * t37) * r_i_i_C(2);
	t42 = r_i_i_C(3) + pkin(9) + pkin(8);
	t39 = t29 * t33;
	t14 = t36 * t26;
	t38 = (-t14 * t21 + t29 * t22) * r_i_i_C(1) + (-t14 * t22 - t29 * t21) * r_i_i_C(2);
	t32 = cos(qJ(4));
	t35 = t32 * pkin(4) + r_i_i_C(1) * t22 - r_i_i_C(2) * t21 + pkin(3);
	t30 = sin(qJ(4));
	t15 = t17 * t29;
	t1 = [0, -t42 * t37 + (-t25 * t39 - t28 * t31) * pkin(2) + t35 * (t25 * t15 - t28 * t36), t41, (t30 * t37 + t32 * t41) * pkin(4) + t43, t43, 0; 0, t42 * t7 + (-t25 * t31 + t28 * t39) * pkin(2) + t35 * (-t28 * t15 - t25 * t36), -t40, (-t30 * t7 - t32 * t40) * pkin(4) + t44, t44, 0; 1, t42 * t14 + (pkin(2) * t33 - t17 * t35) * t26, t29, (-t14 * t30 + t29 * t32) * pkin(4) + t38, t38, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:47
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (366->49), mult. (715->85), div. (0->0), fcn. (934->14), ass. (0->39)
	t67 = pkin(10) + r_i_i_C(3);
	t41 = sin(qJ(6));
	t44 = cos(qJ(6));
	t66 = t44 * r_i_i_C(1) - t41 * r_i_i_C(2) + pkin(5);
	t36 = sin(pkin(12));
	t43 = sin(qJ(2));
	t46 = cos(qJ(2));
	t58 = cos(pkin(12));
	t53 = -t43 * t36 + t46 * t58;
	t35 = qJ(4) + qJ(5);
	t33 = sin(t35);
	t34 = cos(t35);
	t45 = cos(qJ(4));
	t48 = t45 * pkin(4) + t67 * t33 + t66 * t34 + pkin(3);
	t37 = sin(pkin(11));
	t38 = sin(pkin(6));
	t62 = t37 * t38;
	t39 = cos(pkin(11));
	t61 = t39 * t38;
	t40 = cos(pkin(6));
	t60 = t40 * t46;
	t28 = -t46 * t36 - t43 * t58;
	t26 = t28 * t40;
	t13 = -t39 * t26 + t37 * t53;
	t56 = -t37 * t26 - t39 * t53;
	t54 = t41 * r_i_i_C(1) + t44 * r_i_i_C(2) + pkin(8) + pkin(9);
	t8 = t13 * t34 - t33 * t61;
	t52 = t67 * t8 + t66 * (-t13 * t33 - t34 * t61);
	t10 = t33 * t62 - t34 * t56;
	t51 = t67 * t10 + t66 * (t33 * t56 + t34 * t62);
	t50 = t53 * t40;
	t25 = t28 * t38;
	t21 = -t25 * t34 + t40 * t33;
	t49 = t67 * t21 + t66 * (t25 * t33 + t40 * t34);
	t42 = sin(qJ(4));
	t24 = t53 * t38;
	t15 = t39 * t28 - t37 * t50;
	t12 = t37 * t28 + t39 * t50;
	t1 = [0, (-t37 * t60 - t39 * t43) * pkin(2) - t54 * t56 + t48 * t15, t62, (t42 * t56 + t45 * t62) * pkin(4) + t51, t51, (-t10 * t41 - t15 * t44) * r_i_i_C(1) + (-t10 * t44 + t15 * t41) * r_i_i_C(2); 0, (-t37 * t43 + t39 * t60) * pkin(2) + t54 * t13 + t48 * t12, -t61, (-t13 * t42 - t45 * t61) * pkin(4) + t52, t52, (-t12 * t44 - t8 * t41) * r_i_i_C(1) + (t12 * t41 - t8 * t44) * r_i_i_C(2); 1, t38 * t46 * pkin(2) + t48 * t24 - t54 * t25, t40, (t25 * t42 + t40 * t45) * pkin(4) + t49, t49, (-t21 * t41 - t24 * t44) * r_i_i_C(1) + (-t21 * t44 + t24 * t41) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end