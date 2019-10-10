% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:58
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:58:33
	% EndTime: 2019-10-10 12:58:34
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:58:33
	% EndTime: 2019-10-10 12:58:34
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
	% StartTime: 2019-10-10 12:58:33
	% EndTime: 2019-10-10 12:58:34
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
	% StartTime: 2019-10-10 12:58:33
	% EndTime: 2019-10-10 12:58:34
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
	% StartTime: 2019-10-10 12:58:33
	% EndTime: 2019-10-10 12:58:34
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (86->12), mult. (61->18), div. (0->0), fcn. (61->8), ass. (0->15)
	t11 = qJ(2) + qJ(3);
	t8 = qJ(4) + t11;
	t5 = sin(t8);
	t6 = cos(t8);
	t19 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t18 = t19 + pkin(3) * cos(t11);
	t23 = t18 + cos(qJ(2)) * pkin(2);
	t17 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
	t14 = t17 - pkin(3) * sin(t11);
	t20 = r_i_i_C(3) + pkin(9) + pkin(8) + pkin(7);
	t16 = pkin(1) + t23;
	t15 = -sin(qJ(2)) * pkin(2) + t14;
	t13 = cos(qJ(1));
	t12 = sin(qJ(1));
	t1 = [-t16 * t12 + t20 * t13, t15 * t13, t14 * t13, t17 * t13, 0, 0; t20 * t12 + t16 * t13, t15 * t12, t14 * t12, t17 * t12, 0, 0; 0, t23, t18, t19, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:58:33
	% EndTime: 2019-10-10 12:58:34
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (206->33), mult. (162->43), div. (0->0), fcn. (170->10), ass. (0->33)
	t24 = qJ(2) + qJ(3);
	t21 = qJ(4) + t24;
	t18 = sin(t21);
	t19 = cos(t21);
	t25 = sin(qJ(5));
	t44 = r_i_i_C(2) * t25;
	t51 = pkin(10) + r_i_i_C(3);
	t52 = t18 * t44 + t19 * t51;
	t49 = t19 * pkin(4) + t51 * t18;
	t27 = cos(qJ(5));
	t45 = r_i_i_C(1) * t27;
	t33 = (-pkin(4) - t45) * t18;
	t30 = t33 - pkin(3) * sin(t24);
	t17 = pkin(3) * cos(t24);
	t22 = cos(qJ(2)) * pkin(2);
	t48 = pkin(1) + t17 + t22 + t49;
	t26 = sin(qJ(1));
	t41 = t26 * t25;
	t40 = t26 * t27;
	t28 = cos(qJ(1));
	t39 = t28 * t25;
	t38 = t28 * t27;
	t37 = t52 * t26;
	t35 = t52 * t28;
	t32 = -sin(qJ(2)) * pkin(2) + t30;
	t31 = (-t44 + t45) * t19 + t49;
	t29 = t17 + t31;
	t23 = -pkin(9) - pkin(8) - pkin(7);
	t4 = t19 * t38 + t41;
	t3 = -t19 * t39 + t40;
	t2 = -t19 * t40 + t39;
	t1 = t19 * t41 + t38;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t28 * t23 - t26 * t48, t32 * t28 + t35, t30 * t28 + t35, t28 * t33 + t35, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t26 * t23 + t28 * t48, t32 * t26 + t37, t30 * t26 + t37, t26 * t33 + t37, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0; 0, t22 + t29, t29, t31, (-r_i_i_C(1) * t25 - r_i_i_C(2) * t27) * t18, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:58:33
	% EndTime: 2019-10-10 12:58:34
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (306->34), mult. (253->43), div. (0->0), fcn. (274->10), ass. (0->34)
	t28 = sin(qJ(5));
	t30 = cos(qJ(5));
	t39 = r_i_i_C(3) + qJ(6);
	t51 = pkin(5) + r_i_i_C(1);
	t57 = t39 * t28 + t51 * t30;
	t56 = pkin(10) + r_i_i_C(2);
	t27 = qJ(2) + qJ(3);
	t24 = qJ(4) + t27;
	t22 = cos(t24);
	t54 = t22 * t56;
	t21 = sin(t24);
	t53 = t22 * pkin(4) + t56 * t21;
	t34 = (-pkin(4) - t57) * t21;
	t32 = t34 - pkin(3) * sin(t27);
	t18 = pkin(3) * cos(t27);
	t25 = cos(qJ(2)) * pkin(2);
	t52 = pkin(1) + t18 + t25 + t53;
	t29 = sin(qJ(1));
	t49 = t29 * t54;
	t43 = t29 * t28;
	t42 = t29 * t30;
	t31 = cos(qJ(1));
	t41 = t31 * t28;
	t40 = t31 * t30;
	t38 = t31 * t54;
	t36 = t57 * t22 + t53;
	t35 = t18 + t36;
	t33 = -sin(qJ(2)) * pkin(2) + t32;
	t26 = -pkin(9) - pkin(8) - pkin(7);
	t4 = t22 * t40 + t43;
	t3 = t22 * t41 - t42;
	t2 = t22 * t42 - t41;
	t1 = t22 * t43 + t40;
	t5 = [-t39 * t1 - t51 * t2 - t31 * t26 - t52 * t29, t33 * t31 + t38, t32 * t31 + t38, t31 * t34 + t38, -t51 * t3 + t39 * t4, t3; -t29 * t26 + t39 * t3 + t52 * t31 + t51 * t4, t33 * t29 + t49, t32 * t29 + t49, t29 * t34 + t49, -t51 * t1 + t39 * t2, t1; 0, t25 + t35, t35, t36, (-t51 * t28 + t39 * t30) * t21, t21 * t28;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end