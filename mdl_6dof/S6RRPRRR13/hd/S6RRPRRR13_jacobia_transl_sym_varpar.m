% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:09
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR13_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR13_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR13_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:08:56
	% EndTime: 2019-10-10 11:08:56
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:08:56
	% EndTime: 2019-10-10 11:08:56
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
	% StartTime: 2019-10-10 11:08:56
	% EndTime: 2019-10-10 11:08:56
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
	% StartTime: 2019-10-10 11:08:56
	% EndTime: 2019-10-10 11:08:56
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (40->17), mult. (89->25), div. (0->0), fcn. (110->6), ass. (0->18)
	t19 = pkin(2) - r_i_i_C(2);
	t10 = sin(qJ(1));
	t9 = sin(qJ(2));
	t18 = t10 * t9;
	t12 = cos(qJ(1));
	t17 = t12 * t9;
	t11 = cos(qJ(2));
	t16 = t10 * t11;
	t15 = t11 * t12;
	t14 = r_i_i_C(3) + qJ(3);
	t7 = sin(pkin(6));
	t13 = (pkin(8) + r_i_i_C(1)) * t7;
	t8 = cos(pkin(6));
	t4 = -t8 * t18 + t15;
	t3 = t8 * t16 + t17;
	t2 = t8 * t17 + t16;
	t1 = -t8 * t15 + t18;
	t5 = [-t10 * pkin(1) - t14 * t1 + t12 * t13 - t19 * t2, t14 * t4 - t19 * t3, t3, 0, 0, 0; t12 * pkin(1) + t10 * t13 + t14 * t3 + t19 * t4, -t19 * t1 + t14 * t2, t1, 0, 0, 0; 0, (t19 * t11 + t14 * t9) * t7, -t7 * t11, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:08:56
	% EndTime: 2019-10-10 11:08:56
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (75->31), mult. (179->50), div. (0->0), fcn. (222->8), ass. (0->26)
	t27 = pkin(3) + pkin(8);
	t10 = sin(pkin(6));
	t14 = sin(qJ(1));
	t26 = t10 * t14;
	t16 = cos(qJ(2));
	t25 = t10 * t16;
	t17 = cos(qJ(1));
	t24 = t10 * t17;
	t13 = sin(qJ(2));
	t23 = t13 * t14;
	t22 = t13 * t17;
	t21 = t14 * t16;
	t20 = t16 * t17;
	t19 = -r_i_i_C(3) - pkin(9) - pkin(2);
	t12 = sin(qJ(4));
	t15 = cos(qJ(4));
	t18 = r_i_i_C(1) * t12 + r_i_i_C(2) * t15 + qJ(3);
	t11 = cos(pkin(6));
	t8 = t15 * t24;
	t6 = -t11 * t23 + t20;
	t5 = t11 * t21 + t22;
	t4 = t11 * t22 + t21;
	t3 = -t11 * t20 + t23;
	t2 = t12 * t5 + t15 * t26;
	t1 = -t12 * t26 + t15 * t5;
	t7 = [-t14 * pkin(1) + t8 * r_i_i_C(1) - t18 * t3 + (-t12 * r_i_i_C(2) + t27) * t24 + t19 * t4, t18 * t6 + t19 * t5, t5, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0; pkin(1) * t17 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t5 * qJ(3) - t19 * t6 + t27 * t26, t18 * t4 + t19 * t3, t3, (t12 * t24 + t3 * t15) * r_i_i_C(1) + (-t12 * t3 + t8) * r_i_i_C(2), 0, 0; 0, (t18 * t13 - t19 * t16) * t10, -t25, (-t11 * t12 - t15 * t25) * r_i_i_C(1) + (-t11 * t15 + t12 * t25) * r_i_i_C(2), 0, 0;];
	Ja_transl = t7;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:08:56
	% EndTime: 2019-10-10 11:08:56
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (164->47), mult. (410->78), div. (0->0), fcn. (520->10), ass. (0->34)
	t41 = pkin(2) + pkin(9);
	t40 = r_i_i_C(3) + pkin(10);
	t19 = sin(pkin(6));
	t22 = sin(qJ(2));
	t39 = t19 * t22;
	t23 = sin(qJ(1));
	t38 = t19 * t23;
	t26 = cos(qJ(2));
	t37 = t19 * t26;
	t27 = cos(qJ(1));
	t36 = t19 * t27;
	t35 = cos(pkin(6));
	t34 = t19 * (pkin(3) + pkin(8));
	t33 = t23 * t35;
	t32 = t27 * t35;
	t20 = sin(qJ(5));
	t24 = cos(qJ(5));
	t31 = t24 * r_i_i_C(1) - t20 * r_i_i_C(2) + pkin(4);
	t13 = t23 * t22 - t26 * t32;
	t21 = sin(qJ(4));
	t25 = cos(qJ(4));
	t7 = -t13 * t21 + t25 * t36;
	t30 = t13 * t25 + t21 * t36;
	t29 = t20 * r_i_i_C(1) + t24 * r_i_i_C(2) + t41;
	t28 = t31 * t21 - t40 * t25 + qJ(3);
	t16 = -t22 * t33 + t27 * t26;
	t15 = t27 * t22 + t26 * t33;
	t14 = t22 * t32 + t23 * t26;
	t12 = -t21 * t37 + t35 * t25;
	t4 = t15 * t21 + t25 * t38;
	t3 = -t15 * t25 + t21 * t38;
	t2 = t16 * t20 + t4 * t24;
	t1 = t16 * t24 - t4 * t20;
	t5 = [-t23 * pkin(1) - t13 * qJ(3) - t29 * t14 + t27 * t34 + t40 * t30 + t31 * t7, -t29 * t15 + t28 * t16, t15, -t31 * t3 + t40 * t4, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t27 * pkin(1) + t4 * pkin(4) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * qJ(3) + t41 * t16 + t23 * t34 + t40 * t3, -t29 * t13 + t28 * t14, t13, t31 * t30 - t40 * t7, (t14 * t24 + t7 * t20) * r_i_i_C(1) + (-t14 * t20 + t7 * t24) * r_i_i_C(2), 0; 0, (t28 * t22 + t29 * t26) * t19, -t37, t40 * t12 + t31 * (-t35 * t21 - t25 * t37), (-t12 * t20 + t24 * t39) * r_i_i_C(1) + (-t12 * t24 - t20 * t39) * r_i_i_C(2), 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:08:56
	% EndTime: 2019-10-10 11:08:56
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (275->56), mult. (537->89), div. (0->0), fcn. (678->12), ass. (0->41)
	t32 = sin(qJ(2));
	t36 = cos(qJ(2));
	t37 = cos(qJ(1));
	t33 = sin(qJ(1));
	t47 = cos(pkin(6));
	t45 = t33 * t47;
	t21 = t37 * t32 + t36 * t45;
	t31 = sin(qJ(4));
	t35 = cos(qJ(4));
	t29 = sin(pkin(6));
	t50 = t29 * t33;
	t10 = t21 * t31 + t35 * t50;
	t22 = -t32 * t45 + t37 * t36;
	t28 = qJ(5) + qJ(6);
	t26 = sin(t28);
	t27 = cos(t28);
	t5 = -t10 * t26 + t22 * t27;
	t6 = t10 * t27 + t22 * t26;
	t55 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t44 = t37 * t47;
	t19 = t33 * t32 - t36 * t44;
	t48 = t29 * t37;
	t13 = -t19 * t31 + t35 * t48;
	t20 = t32 * t44 + t33 * t36;
	t54 = (t13 * t26 + t20 * t27) * r_i_i_C(1) + (t13 * t27 - t20 * t26) * r_i_i_C(2);
	t49 = t29 * t36;
	t18 = -t31 * t49 + t47 * t35;
	t51 = t29 * t32;
	t53 = (-t18 * t26 + t27 * t51) * r_i_i_C(1) + (-t18 * t27 - t26 * t51) * r_i_i_C(2);
	t52 = r_i_i_C(3) + pkin(11) + pkin(10);
	t46 = t29 * (pkin(3) + pkin(8));
	t30 = sin(qJ(5));
	t43 = t30 * pkin(5) + pkin(2) + pkin(9);
	t34 = cos(qJ(5));
	t25 = t34 * pkin(5) + pkin(4);
	t42 = r_i_i_C(1) * t27 - r_i_i_C(2) * t26 + t25;
	t41 = t19 * t35 + t31 * t48;
	t40 = t26 * r_i_i_C(1) + t27 * r_i_i_C(2) + t43;
	t39 = t42 * t31 - t52 * t35 + qJ(3);
	t9 = -t21 * t35 + t31 * t50;
	t1 = [-t33 * pkin(1) - t19 * qJ(3) + t42 * t13 - t40 * t20 + t37 * t46 + t52 * t41, -t40 * t21 + t39 * t22, t21, t52 * t10 - t42 * t9, (-t10 * t30 + t22 * t34) * pkin(5) + t55, t55; t37 * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t21 * qJ(3) + t10 * t25 + t43 * t22 + t33 * t46 + t52 * t9, -t40 * t19 + t39 * t20, t19, -t52 * t13 + t42 * t41, (t13 * t30 + t20 * t34) * pkin(5) + t54, t54; 0, (t39 * t32 + t40 * t36) * t29, -t49, t52 * t18 + t42 * (-t47 * t31 - t35 * t49), (-t18 * t30 + t34 * t51) * pkin(5) + t53, t53;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end