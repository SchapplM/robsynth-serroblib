% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR12
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
% Datum: 2019-10-10 11:07
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR12_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR12_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR12_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:01
	% EndTime: 2019-10-10 11:07:01
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:01
	% EndTime: 2019-10-10 11:07:01
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
	% StartTime: 2019-10-10 11:07:01
	% EndTime: 2019-10-10 11:07:01
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
	% StartTime: 2019-10-10 11:07:01
	% EndTime: 2019-10-10 11:07:01
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
	% StartTime: 2019-10-10 11:07:01
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.15s
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
	t18 = t12 * r_i_i_C(1) + t15 * r_i_i_C(2) + qJ(3);
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
	% StartTime: 2019-10-10 11:07:01
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (150->40), mult. (261->61), div. (0->0), fcn. (320->10), ass. (0->33)
	t21 = cos(pkin(6));
	t23 = sin(qJ(2));
	t27 = cos(qJ(1));
	t33 = t27 * t23;
	t24 = sin(qJ(1));
	t26 = cos(qJ(2));
	t34 = t24 * t26;
	t11 = t21 * t34 + t33;
	t19 = qJ(4) + qJ(5);
	t17 = sin(t19);
	t18 = cos(t19);
	t20 = sin(pkin(6));
	t38 = t20 * t24;
	t5 = t11 * t18 - t17 * t38;
	t6 = t11 * t17 + t18 * t38;
	t42 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t36 = t20 * t27;
	t13 = t18 * t36;
	t32 = t27 * t26;
	t35 = t24 * t23;
	t9 = -t21 * t32 + t35;
	t41 = (t17 * t36 + t9 * t18) * r_i_i_C(1) + (-t9 * t17 + t13) * r_i_i_C(2);
	t37 = t20 * t26;
	t40 = (-t21 * t17 - t18 * t37) * r_i_i_C(1) + (t17 * t37 - t21 * t18) * r_i_i_C(2);
	t25 = cos(qJ(4));
	t39 = t25 * pkin(4) + pkin(3) + pkin(8);
	t31 = -r_i_i_C(3) - pkin(10) - pkin(9) - pkin(2);
	t22 = sin(qJ(4));
	t30 = t22 * pkin(4) + qJ(3);
	t29 = t17 * r_i_i_C(1) + t18 * r_i_i_C(2) + t30;
	t12 = -t21 * t35 + t32;
	t10 = t21 * t33 + t34;
	t1 = [-t24 * pkin(1) + t13 * r_i_i_C(1) - t29 * t9 + (-r_i_i_C(2) * t17 + t39) * t36 + t31 * t10, t11 * t31 + t12 * t29, t11, (t11 * t25 - t22 * t38) * pkin(4) + t42, t42, 0; t27 * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t11 * t30 - t12 * t31 + t38 * t39, t10 * t29 + t31 * t9, t9, (t22 * t36 + t25 * t9) * pkin(4) + t41, t41, 0; 0, (t23 * t29 - t26 * t31) * t20, -t37, (-t21 * t22 - t25 * t37) * pkin(4) + t40, t40, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:01
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (332->56), mult. (552->89), div. (0->0), fcn. (693->12), ass. (0->41)
	t61 = r_i_i_C(3) + pkin(11);
	t34 = sin(qJ(6));
	t38 = cos(qJ(6));
	t48 = t38 * r_i_i_C(1) - t34 * r_i_i_C(2) + pkin(5);
	t58 = pkin(2) + pkin(10) + pkin(9);
	t33 = sin(pkin(6));
	t36 = sin(qJ(2));
	t57 = t33 * t36;
	t37 = sin(qJ(1));
	t56 = t33 * t37;
	t40 = cos(qJ(2));
	t55 = t33 * t40;
	t41 = cos(qJ(1));
	t54 = t33 * t41;
	t53 = cos(pkin(6));
	t39 = cos(qJ(4));
	t52 = t33 * (t39 * pkin(4) + pkin(3) + pkin(8));
	t35 = sin(qJ(4));
	t51 = t35 * pkin(4) + qJ(3);
	t50 = t37 * t53;
	t49 = t41 * t53;
	t22 = t37 * t36 - t40 * t49;
	t32 = qJ(4) + qJ(5);
	t30 = sin(t32);
	t31 = cos(t32);
	t13 = -t22 * t30 + t31 * t54;
	t11 = t22 * t31 + t30 * t54;
	t47 = t34 * r_i_i_C(1) + t38 * r_i_i_C(2) + t58;
	t24 = t41 * t36 + t40 * t50;
	t10 = t24 * t30 + t31 * t56;
	t9 = -t24 * t31 + t30 * t56;
	t46 = t61 * t10 - t48 * t9;
	t45 = t48 * t11 - t61 * t13;
	t21 = -t30 * t55 + t53 * t31;
	t44 = t61 * t21 + t48 * (-t53 * t30 - t31 * t55);
	t43 = t48 * t30 - t61 * t31 + t51;
	t25 = -t36 * t50 + t41 * t40;
	t23 = t36 * t49 + t37 * t40;
	t2 = t10 * t38 + t25 * t34;
	t1 = -t10 * t34 + t25 * t38;
	t3 = [-t37 * pkin(1) + t61 * t11 + t48 * t13 - t51 * t22 - t47 * t23 + t41 * t52, -t47 * t24 + t43 * t25, t24, (t24 * t39 - t35 * t56) * pkin(4) + t46, t46, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t41 * pkin(1) + t10 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t51 * t24 + t58 * t25 + t37 * t52 + t61 * t9, -t47 * t22 + t43 * t23, t22, (t22 * t39 + t35 * t54) * pkin(4) + t45, t45, (t13 * t34 + t23 * t38) * r_i_i_C(1) + (t13 * t38 - t23 * t34) * r_i_i_C(2); 0, (t43 * t36 + t47 * t40) * t33, -t55, (-t53 * t35 - t39 * t55) * pkin(4) + t44, t44, (-t21 * t34 + t38 * t57) * r_i_i_C(1) + (-t21 * t38 - t34 * t57) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end