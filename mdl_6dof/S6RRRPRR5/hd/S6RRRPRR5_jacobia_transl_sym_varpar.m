% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:00
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
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
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
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
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
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
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (73->17), mult. (71->20), div. (0->0), fcn. (74->6), ass. (0->17)
	t33 = r_i_i_C(3) + qJ(4);
	t14 = qJ(2) + qJ(3);
	t11 = sin(t14);
	t12 = cos(t14);
	t19 = (pkin(3) - r_i_i_C(2)) * t12 + t33 * t11;
	t32 = cos(qJ(2)) * pkin(2) + t19;
	t31 = t33 * t12;
	t30 = pkin(1) + t32;
	t27 = r_i_i_C(1) + pkin(8) + pkin(7);
	t16 = sin(qJ(1));
	t26 = t16 * t11;
	t17 = cos(qJ(1));
	t25 = t17 * t11;
	t22 = r_i_i_C(2) * t26 + t31 * t16;
	t21 = r_i_i_C(2) * t25 + t31 * t17;
	t20 = -sin(qJ(2)) * pkin(2) - pkin(3) * t11;
	t1 = [-t30 * t16 + t27 * t17, t20 * t17 + t21, -pkin(3) * t25 + t21, t25, 0, 0; t27 * t16 + t30 * t17, t20 * t16 + t22, -pkin(3) * t26 + t22, t26, 0, 0; 0, t32, t19, -t12, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (122->28), mult. (139->40), div. (0->0), fcn. (150->8), ass. (0->27)
	t22 = sin(qJ(5));
	t25 = cos(qJ(5));
	t48 = r_i_i_C(1) * t22 + r_i_i_C(2) * t25;
	t21 = qJ(2) + qJ(3);
	t18 = sin(t21);
	t19 = cos(t21);
	t36 = pkin(3) + pkin(9) + r_i_i_C(3);
	t47 = t18 * qJ(4) + t36 * t19;
	t45 = (qJ(4) + t48) * t19;
	t20 = cos(qJ(2)) * pkin(2);
	t44 = t20 + pkin(1) + t47;
	t41 = pkin(4) + pkin(8) + pkin(7);
	t24 = sin(qJ(1));
	t40 = t24 * t18;
	t39 = t24 * t25;
	t26 = cos(qJ(1));
	t38 = t26 * t18;
	t35 = t45 * t24;
	t34 = t45 * t26;
	t30 = t36 * t18;
	t29 = t48 * t18 + t47;
	t28 = -sin(qJ(2)) * pkin(2) - t30;
	t4 = -t22 * t40 + t25 * t26;
	t3 = t18 * t39 + t22 * t26;
	t2 = t22 * t38 + t39;
	t1 = -t22 * t24 + t25 * t38;
	t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t44 * t24 + t41 * t26, t28 * t26 + t34, -t26 * t30 + t34, t38, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t41 * t24 + t44 * t26, t28 * t24 + t35, -t24 * t30 + t35, t40, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0; 0, t20 + t29, t29, -t19, (-r_i_i_C(1) * t25 + r_i_i_C(2) * t22) * t19, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (209->43), mult. (195->55), div. (0->0), fcn. (210->10), ass. (0->38)
	t35 = sin(qJ(5));
	t67 = pkin(5) * t35 + qJ(4);
	t34 = qJ(2) + qJ(3);
	t31 = cos(t34);
	t40 = -pkin(10) - pkin(9);
	t64 = -pkin(3) - r_i_i_C(3);
	t68 = (-t40 - t64) * t31;
	t33 = qJ(5) + qJ(6);
	t30 = cos(t33);
	t58 = t30 * t31;
	t28 = sin(t33);
	t59 = t28 * t31;
	t66 = r_i_i_C(1) * t59 + r_i_i_C(2) * t58 + t67 * t31;
	t29 = sin(t34);
	t32 = cos(qJ(2)) * pkin(2);
	t65 = t67 * t29 + pkin(1) + t32 + t68;
	t37 = sin(qJ(1));
	t39 = cos(qJ(1));
	t54 = t39 * t30;
	t5 = -t37 * t28 + t29 * t54;
	t55 = t39 * t29;
	t56 = t37 * t30;
	t6 = t28 * t55 + t56;
	t63 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t7 = t39 * t28 + t29 * t56;
	t57 = t37 * t29;
	t8 = -t28 * t57 + t54;
	t62 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
	t38 = cos(qJ(5));
	t60 = t38 * pkin(5);
	t53 = pkin(4) + t60 + pkin(8) + pkin(7);
	t47 = t64 * t29;
	t45 = t66 * t37 + t40 * t57;
	t44 = t66 * t39 + t40 * t55;
	t43 = -sin(qJ(2)) * pkin(2) + t47;
	t42 = t68 + (r_i_i_C(1) * t28 + r_i_i_C(2) * t30 + t67) * t29;
	t17 = r_i_i_C(2) * t59;
	t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) - t65 * t37 + t53 * t39, t43 * t39 + t44, t39 * t47 + t44, t55, (-t35 * t37 + t38 * t55) * pkin(5) + t63, t63; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t53 * t37 + t65 * t39, t43 * t37 + t45, t37 * t47 + t45, t57, (t35 * t39 + t38 * t57) * pkin(5) + t62, t62; 0, t32 + t42, t42, -t31, t17 + (-r_i_i_C(1) * t30 - t60) * t31, -r_i_i_C(1) * t58 + t17;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end