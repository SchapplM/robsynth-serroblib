% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR14V3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14V3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR14V3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobia_transl_sym_varpar: pkin has to be [1x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:22
	% EndTime: 2019-10-10 11:13:22
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:22
	% EndTime: 2019-10-10 11:13:23
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
	% StartTime: 2019-10-10 11:13:22
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (7->4), mult. (20->10), div. (0->0), fcn. (20->4), ass. (0->7)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t6 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t5 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t7 = [t4 * r_i_i_C(3) - t2 * t6, t5 * t4, 0, 0, 0, 0; t2 * r_i_i_C(3) + t4 * t6, t5 * t2, 0, 0, 0, 0; 0, t6, 0, 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:22
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (13->6), mult. (31->12), div. (0->0), fcn. (34->4), ass. (0->8)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(3) + qJ(3);
	t5 = t3 * r_i_i_C(1) + t7 * t1;
	t6 = -r_i_i_C(1) * t1 + t7 * t3;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t8 = [t4 * r_i_i_C(2) - t5 * t2, t6 * t4, t4 * t1, 0, 0, 0; t2 * r_i_i_C(2) + t5 * t4, t6 * t2, t2 * t1, 0, 0, 0; 0, t5, -t3, 0, 0, 0;];
	Ja_transl = t8;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:22
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (27->16), mult. (72->33), div. (0->0), fcn. (83->6), ass. (0->18)
	t6 = sin(qJ(2));
	t7 = sin(qJ(1));
	t17 = t7 * t6;
	t9 = cos(qJ(2));
	t16 = t7 * t9;
	t10 = cos(qJ(1));
	t15 = t10 * t9;
	t14 = r_i_i_C(3) + qJ(3);
	t13 = t14 * t6;
	t5 = sin(qJ(4));
	t8 = cos(qJ(4));
	t12 = r_i_i_C(1) * t8 - r_i_i_C(2) * t5;
	t11 = -t12 * t6 + t14 * t9;
	t4 = t8 * t15 + t7 * t5;
	t3 = -t5 * t15 + t7 * t8;
	t2 = t10 * t5 - t8 * t16;
	t1 = t10 * t8 + t5 * t16;
	t18 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t14 * t17, t11 * t10, t10 * t6, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t10 * t13, t11 * t7, t17, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0; 0, t12 * t9 + t13, -t9, (-r_i_i_C(1) * t5 - r_i_i_C(2) * t8) * t6, 0, 0;];
	Ja_transl = t18;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:22
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (57->32), mult. (159->62), div. (0->0), fcn. (194->8), ass. (0->30)
	t10 = sin(qJ(4));
	t31 = t10 * r_i_i_C(3);
	t11 = sin(qJ(2));
	t14 = cos(qJ(4));
	t30 = t11 * t14;
	t12 = sin(qJ(1));
	t29 = t12 * t10;
	t28 = t12 * t11;
	t15 = cos(qJ(2));
	t27 = t14 * t15;
	t16 = cos(qJ(1));
	t26 = t14 * t16;
	t25 = t16 * t10;
	t24 = t16 * t11;
	t23 = t11 * qJ(3);
	t13 = cos(qJ(5));
	t9 = sin(qJ(5));
	t22 = t13 * r_i_i_C(1) - t9 * r_i_i_C(2);
	t4 = t12 * t27 - t25;
	t21 = -t13 * t4 - t9 * t28;
	t20 = t13 * t28 - t4 * t9;
	t19 = t13 * t15 + t9 * t30;
	t18 = -t13 * t30 + t15 * t9;
	t17 = t18 * r_i_i_C(1) + t19 * r_i_i_C(2) + t15 * qJ(3) - t11 * t31;
	t6 = t15 * t26 + t29;
	t5 = -t12 * t14 + t15 * t25;
	t3 = -t15 * t29 - t26;
	t2 = t6 * t13 + t9 * t24;
	t1 = t13 * t24 - t6 * t9;
	t7 = [t21 * r_i_i_C(1) - t20 * r_i_i_C(2) + t3 * r_i_i_C(3) - t12 * t23, t17 * t16, t24, r_i_i_C(3) * t6 - t22 * t5, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t5 * r_i_i_C(3) + t16 * t23, t17 * t12, t28, r_i_i_C(3) * t4 + t22 * t3, t20 * r_i_i_C(1) + t21 * r_i_i_C(2), 0; 0, (t11 * t9 + t13 * t27) * r_i_i_C(1) + (t11 * t13 - t9 * t27) * r_i_i_C(2) + t15 * t31 + t23, -t15, (r_i_i_C(3) * t14 - t22 * t10) * t11, -t19 * r_i_i_C(1) + t18 * r_i_i_C(2), 0;];
	Ja_transl = t7;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:22
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (117->53), mult. (330->105), div. (0->0), fcn. (419->10), ass. (0->42)
	t23 = sin(qJ(4));
	t28 = cos(qJ(4));
	t30 = cos(qJ(1));
	t37 = t30 * t28;
	t25 = sin(qJ(1));
	t29 = cos(qJ(2));
	t43 = t25 * t29;
	t12 = t23 * t43 + t37;
	t21 = sin(qJ(6));
	t26 = cos(qJ(6));
	t39 = t30 * t23;
	t13 = t28 * t43 - t39;
	t22 = sin(qJ(5));
	t27 = cos(qJ(5));
	t24 = sin(qJ(2));
	t44 = t25 * t24;
	t4 = t13 * t27 + t22 * t44;
	t54 = -t12 * t26 + t4 * t21;
	t53 = -t12 * t21 - t4 * t26;
	t41 = t29 * t22;
	t45 = t24 * t28;
	t11 = t27 * t45 - t41;
	t40 = t29 * t27;
	t33 = -t22 * t45 - t40;
	t34 = t26 * r_i_i_C(1) - t21 * r_i_i_C(2);
	t47 = t23 * t24;
	t52 = t33 * r_i_i_C(3) + t29 * qJ(3) + (-t21 * r_i_i_C(1) - t26 * r_i_i_C(2)) * t47 - t34 * t11;
	t51 = t22 * r_i_i_C(3);
	t48 = t21 * t27;
	t46 = t23 * t29;
	t42 = t26 * t27;
	t38 = t30 * t24;
	t36 = t24 * qJ(3);
	t35 = -t13 * t22 + t27 * t44;
	t16 = t25 * t23 + t29 * t37;
	t15 = -t25 * t28 + t29 * t39;
	t14 = t24 * t22 + t28 * t40;
	t7 = t16 * t27 + t22 * t38;
	t6 = t16 * t22 - t27 * t38;
	t2 = t15 * t21 + t7 * t26;
	t1 = t15 * t26 - t7 * t21;
	t3 = [t53 * r_i_i_C(1) + t54 * r_i_i_C(2) + t35 * r_i_i_C(3) - t25 * t36, t52 * t30, t38, (-t15 * t42 + t16 * t21) * r_i_i_C(1) + (t15 * t48 + t16 * t26) * r_i_i_C(2) - t15 * t51, t7 * r_i_i_C(3) - t34 * t6, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t6 * r_i_i_C(3) + t30 * t36, t52 * t25, t44, (-t12 * t42 + t13 * t21) * r_i_i_C(1) + (t12 * t48 + t13 * t26) * r_i_i_C(2) - t12 * t51, t4 * r_i_i_C(3) + t34 * t35, -t54 * r_i_i_C(1) + t53 * r_i_i_C(2); 0, (t14 * t26 + t21 * t46) * r_i_i_C(1) + (-t14 * t21 + t26 * t46) * r_i_i_C(2) + (-t24 * t27 + t28 * t41) * r_i_i_C(3) + t36, -t29, ((t21 * t28 - t23 * t42) * r_i_i_C(1) + (t23 * t48 + t26 * t28) * r_i_i_C(2) - t23 * t51) * t24, t11 * r_i_i_C(3) + t34 * t33, (-t11 * t21 + t26 * t47) * r_i_i_C(1) + (-t11 * t26 - t21 * t47) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end