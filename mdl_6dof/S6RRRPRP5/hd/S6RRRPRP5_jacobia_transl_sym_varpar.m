% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:42
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:42:10
	% EndTime: 2019-10-10 11:42:11
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:42:10
	% EndTime: 2019-10-10 11:42:11
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
	% StartTime: 2019-10-10 11:42:11
	% EndTime: 2019-10-10 11:42:11
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
	% StartTime: 2019-10-10 11:42:10
	% EndTime: 2019-10-10 11:42:11
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (35->20), mult. (83->34), div. (0->0), fcn. (91->6), ass. (0->18)
	t16 = pkin(8) + r_i_i_C(3);
	t6 = sin(qJ(2));
	t14 = t16 * t6;
	t9 = cos(qJ(2));
	t18 = t9 * pkin(2) + pkin(1) + t14;
	t7 = sin(qJ(1));
	t17 = t7 * t9;
	t10 = cos(qJ(1));
	t15 = t10 * t9;
	t5 = sin(qJ(3));
	t8 = cos(qJ(3));
	t12 = r_i_i_C(1) * t8 - r_i_i_C(2) * t5 + pkin(2);
	t11 = -t12 * t6 + t16 * t9;
	t4 = t8 * t15 + t7 * t5;
	t3 = -t5 * t15 + t7 * t8;
	t2 = t10 * t5 - t8 * t17;
	t1 = t10 * t8 + t5 * t17;
	t13 = [t10 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t18 * t7, t11 * t10, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0, 0; t7 * pkin(7) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t18 * t10, t11 * t7, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0, 0; 0, t12 * t9 + t14, (-r_i_i_C(1) * t5 - r_i_i_C(2) * t8) * t6, 0, 0, 0;];
	Ja_transl = t13;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:42:10
	% EndTime: 2019-10-10 11:42:11
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (77->29), mult. (106->44), div. (0->0), fcn. (117->8), ass. (0->24)
	t14 = cos(qJ(2));
	t11 = sin(qJ(2));
	t24 = r_i_i_C(3) + qJ(4) + pkin(8);
	t20 = t24 * t11;
	t13 = cos(qJ(3));
	t5 = pkin(3) * t13 + pkin(2);
	t25 = t14 * t5 + pkin(1) + t20;
	t10 = sin(qJ(3));
	t23 = pkin(3) * t10;
	t12 = sin(qJ(1));
	t22 = t12 * t14;
	t15 = cos(qJ(1));
	t21 = t14 * t15;
	t18 = pkin(7) + t23;
	t8 = qJ(3) + pkin(10);
	t6 = sin(t8);
	t7 = cos(t8);
	t17 = r_i_i_C(1) * t7 - r_i_i_C(2) * t6 + t5;
	t16 = -t17 * t11 + t24 * t14;
	t4 = t12 * t6 + t7 * t21;
	t3 = t12 * t7 - t6 * t21;
	t2 = t15 * t6 - t7 * t22;
	t1 = t15 * t7 + t6 * t22;
	t9 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t25 * t12 + t18 * t15, t16 * t15, t3 * r_i_i_C(1) - t4 * r_i_i_C(2) + (-t10 * t21 + t12 * t13) * pkin(3), t15 * t11, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t18 * t12 + t25 * t15, t16 * t12, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + (-t10 * t22 - t13 * t15) * pkin(3), t12 * t11, 0, 0; 0, t17 * t14 + t20, (-r_i_i_C(1) * t6 - r_i_i_C(2) * t7 - t23) * t11, -t14, 0, 0;];
	Ja_transl = t9;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:42:10
	% EndTime: 2019-10-10 11:42:11
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (157->32), mult. (138->45), div. (0->0), fcn. (153->10), ass. (0->28)
	t21 = cos(qJ(2));
	t19 = sin(qJ(2));
	t31 = r_i_i_C(3) + pkin(9) + qJ(4) + pkin(8);
	t26 = t31 * t19;
	t18 = qJ(3) + pkin(10);
	t11 = pkin(4) * cos(t18) + cos(qJ(3)) * pkin(3);
	t9 = pkin(2) + t11;
	t35 = t21 * t9 + pkin(1) + t26;
	t15 = qJ(5) + t18;
	t13 = sin(t15);
	t14 = cos(t15);
	t22 = cos(qJ(1));
	t28 = t22 * t14;
	t20 = sin(qJ(1));
	t30 = t20 * t21;
	t5 = t13 * t30 + t28;
	t29 = t22 * t13;
	t6 = -t14 * t30 + t29;
	t34 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t7 = t20 * t14 - t21 * t29;
	t8 = t20 * t13 + t21 * t28;
	t33 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t10 = pkin(4) * sin(t18) + sin(qJ(3)) * pkin(3);
	t32 = pkin(7) + t10;
	t25 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
	t24 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t9;
	t23 = -t24 * t19 + t31 * t21;
	t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t35 * t20 + t32 * t22, t23 * t22, -t22 * t21 * t10 + t20 * t11 + t33, t22 * t19, t33, 0; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t32 * t20 + t35 * t22, t23 * t20, -t10 * t30 - t22 * t11 + t34, t20 * t19, t34, 0; 0, t24 * t21 + t26, (-t10 + t25) * t19, -t21, t25 * t19, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:42:10
	% EndTime: 2019-10-10 11:42:11
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (271->35), mult. (217->46), div. (0->0), fcn. (249->10), ass. (0->31)
	t43 = pkin(5) + r_i_i_C(1);
	t36 = r_i_i_C(3) + qJ(6);
	t24 = qJ(3) + pkin(10);
	t13 = pkin(4) * cos(t24) + cos(qJ(3)) * pkin(3);
	t11 = pkin(2) + t13;
	t27 = cos(qJ(2));
	t25 = sin(qJ(2));
	t41 = r_i_i_C(2) + pkin(9) + qJ(4) + pkin(8);
	t33 = t41 * t25;
	t45 = t27 * t11 + pkin(1) + t33;
	t21 = qJ(5) + t24;
	t19 = sin(t21);
	t20 = cos(t21);
	t44 = t36 * t19 + t43 * t20 + t11;
	t12 = pkin(4) * sin(t24) + sin(qJ(3)) * pkin(3);
	t42 = pkin(7) + t12;
	t26 = sin(qJ(1));
	t39 = t26 * t27;
	t28 = cos(qJ(1));
	t38 = t28 * t19;
	t37 = t28 * t20;
	t35 = t36 * t20 * t25;
	t34 = t43 * t19;
	t7 = t19 * t39 + t37;
	t8 = t20 * t39 - t38;
	t31 = t36 * t8 - t43 * t7;
	t10 = t26 * t19 + t27 * t37;
	t9 = -t26 * t20 + t27 * t38;
	t30 = t36 * t10 - t43 * t9;
	t29 = -t44 * t25 + t41 * t27;
	t1 = [-t45 * t26 + t42 * t28 - t36 * t7 - t43 * t8, t29 * t28, -t28 * t27 * t12 + t26 * t13 + t30, t28 * t25, t30, t9; t43 * t10 + t42 * t26 + t45 * t28 + t36 * t9, t29 * t26, -t12 * t39 - t28 * t13 + t31, t26 * t25, t31, t7; 0, t44 * t27 + t33, (-t12 - t34) * t25 + t35, -t27, -t25 * t34 + t35, t25 * t19;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end