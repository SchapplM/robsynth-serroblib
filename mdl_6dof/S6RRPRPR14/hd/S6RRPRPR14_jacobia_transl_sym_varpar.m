% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:28
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR14_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR14_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR14_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:15
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
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:15
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
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:15
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
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:15
	% DurationCPUTime: 0.21s
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
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:15
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (118->34), mult. (284->51), div. (0->0), fcn. (357->8), ass. (0->30)
	t36 = pkin(4) - r_i_i_C(2);
	t16 = sin(pkin(6));
	t20 = sin(qJ(1));
	t35 = t16 * t20;
	t22 = cos(qJ(2));
	t34 = t16 * t22;
	t23 = cos(qJ(1));
	t33 = t16 * t23;
	t19 = sin(qJ(2));
	t32 = t20 * t19;
	t31 = t20 * t22;
	t30 = t23 * t19;
	t29 = t23 * t22;
	t28 = r_i_i_C(3) + qJ(5);
	t27 = pkin(2) + pkin(9) + r_i_i_C(1);
	t26 = t16 * (pkin(3) + pkin(8));
	t17 = cos(pkin(6));
	t10 = -t17 * t29 + t32;
	t18 = sin(qJ(4));
	t21 = cos(qJ(4));
	t25 = -t10 * t18 + t21 * t33;
	t3 = t10 * t21 + t18 * t33;
	t24 = t36 * t18 - t28 * t21 + qJ(3);
	t13 = -t17 * t32 + t29;
	t12 = t17 * t31 + t30;
	t11 = t17 * t30 + t31;
	t8 = t17 * t18 + t21 * t34;
	t2 = t12 * t18 + t21 * t35;
	t1 = -t12 * t21 + t18 * t35;
	t4 = [-t20 * pkin(1) - t10 * qJ(3) - t27 * t11 + t23 * t26 + t36 * t25 + t28 * t3, -t27 * t12 + t24 * t13, t12, -t36 * t1 + t28 * t2, t1, 0; t23 * pkin(1) + t12 * qJ(3) + t28 * t1 + t27 * t13 + t36 * t2 + t20 * t26, -t27 * t10 + t24 * t11, t10, -t28 * t25 + t36 * t3, -t3, 0; 0, (t24 * t19 + t27 * t22) * t16, -t34, t28 * (t17 * t21 - t18 * t34) - t36 * t8, t8, 0;];
	Ja_transl = t4;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:15
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (197->48), mult. (488->78), div. (0->0), fcn. (621->10), ass. (0->34)
	t19 = sin(pkin(6));
	t22 = sin(qJ(2));
	t41 = t19 * t22;
	t23 = sin(qJ(1));
	t40 = t19 * t23;
	t26 = cos(qJ(2));
	t39 = t19 * t26;
	t27 = cos(qJ(1));
	t38 = t19 * t27;
	t37 = cos(pkin(6));
	t36 = pkin(2) + pkin(5) + pkin(9);
	t35 = pkin(4) + pkin(10) + r_i_i_C(3);
	t34 = t19 * (pkin(3) + pkin(8));
	t33 = t23 * t37;
	t32 = t27 * t37;
	t20 = sin(qJ(6));
	t24 = cos(qJ(6));
	t31 = t20 * r_i_i_C(1) + t24 * r_i_i_C(2) + qJ(5);
	t13 = t23 * t22 - t26 * t32;
	t21 = sin(qJ(4));
	t25 = cos(qJ(4));
	t30 = -t13 * t21 + t25 * t38;
	t5 = t13 * t25 + t21 * t38;
	t29 = t24 * r_i_i_C(1) - t20 * r_i_i_C(2) + t36;
	t28 = t35 * t21 - t31 * t25 + qJ(3);
	t16 = -t22 * t33 + t27 * t26;
	t15 = t27 * t22 + t26 * t33;
	t14 = t22 * t32 + t23 * t26;
	t11 = t37 * t21 + t25 * t39;
	t4 = t15 * t21 + t25 * t40;
	t3 = -t15 * t25 + t21 * t40;
	t2 = t16 * t24 + t3 * t20;
	t1 = -t16 * t20 + t3 * t24;
	t6 = [-t23 * pkin(1) - t13 * qJ(3) - t29 * t14 + t27 * t34 + t35 * t30 + t31 * t5, -t29 * t15 + t28 * t16, t15, -t35 * t3 + t31 * t4, t3, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t27 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * qJ(3) + t3 * qJ(5) + t36 * t16 + t23 * t34 + t35 * t4, -t29 * t13 + t28 * t14, t13, -t31 * t30 + t35 * t5, -t5, (-t14 * t20 - t5 * t24) * r_i_i_C(1) + (-t14 * t24 + t5 * t20) * r_i_i_C(2); 0, (t28 * t22 + t29 * t26) * t19, -t39, t31 * (-t21 * t39 + t37 * t25) - t35 * t11, t11, (t11 * t24 - t20 * t41) * r_i_i_C(1) + (-t11 * t20 - t24 * t41) * r_i_i_C(2);];
	Ja_transl = t6;
else
	Ja_transl=NaN(3,6);
end