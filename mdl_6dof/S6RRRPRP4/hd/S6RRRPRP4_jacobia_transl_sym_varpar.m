% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:40
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:25
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:25
	% EndTime: 2019-10-10 11:40:26
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
	% StartTime: 2019-10-10 11:40:25
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (11->6), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->9)
	t8 = pkin(7) + r_i_i_C(3);
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = t3 * r_i_i_C(1) - t1 * r_i_i_C(2);
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t5 * t2 + t8 * t4, t6 * t4, 0, 0, 0, 0; t8 * t2 + t5 * t4, t6 * t2, 0, 0, 0, 0; 0, t7, 0, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:25
	% EndTime: 2019-10-10 11:40:26
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
	% StartTime: 2019-10-10 11:40:25
	% EndTime: 2019-10-10 11:40:26
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
	% StartTime: 2019-10-10 11:40:25
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.22s
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
	% StartTime: 2019-10-10 11:40:25
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (179->33), mult. (214->45), div. (0->0), fcn. (238->8), ass. (0->31)
	t22 = qJ(2) + qJ(3);
	t19 = sin(t22);
	t20 = cos(t22);
	t38 = pkin(3) + pkin(9) + r_i_i_C(2);
	t51 = t19 * qJ(4) + t38 * t20;
	t47 = pkin(5) + r_i_i_C(1);
	t49 = t20 * t47;
	t21 = cos(qJ(2)) * pkin(2);
	t48 = t21 + pkin(1) + t51;
	t46 = pkin(4) + pkin(8) + pkin(7);
	t23 = sin(qJ(5));
	t25 = sin(qJ(1));
	t44 = t25 * t23;
	t26 = cos(qJ(5));
	t43 = t25 * t26;
	t27 = cos(qJ(1));
	t42 = t27 * t23;
	t41 = t27 * t26;
	t40 = r_i_i_C(3) + qJ(6);
	t39 = qJ(4) * t20;
	t37 = t25 * t39 + t44 * t49;
	t36 = t27 * t39 + t42 * t49;
	t33 = t40 * t26;
	t31 = -t38 * t19 - t20 * t33;
	t30 = (t23 * t47 - t33) * t19 + t51;
	t29 = -sin(qJ(2)) * pkin(2) + t31;
	t4 = -t19 * t44 + t41;
	t3 = t19 * t43 + t42;
	t2 = t19 * t42 + t43;
	t1 = -t19 * t41 + t44;
	t5 = [-t48 * t25 + t46 * t27 + t40 * t3 + t47 * t4, t29 * t27 + t36, t31 * t27 + t36, t27 * t19, -t47 * t1 + t40 * t2, t1; t40 * t1 + t47 * t2 + t46 * t25 + t48 * t27, t29 * t25 + t37, t31 * t25 + t37, t25 * t19, t47 * t3 - t40 * t4, -t3; 0, t21 + t30, t30, -t20, (-t40 * t23 - t47 * t26) * t20, t20 * t26;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end