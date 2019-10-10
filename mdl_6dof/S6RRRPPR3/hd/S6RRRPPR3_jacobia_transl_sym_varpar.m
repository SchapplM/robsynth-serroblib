% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:20
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:20:30
	% EndTime: 2019-10-10 11:20:30
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:20:30
	% EndTime: 2019-10-10 11:20:30
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
	% StartTime: 2019-10-10 11:20:30
	% EndTime: 2019-10-10 11:20:30
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
	% StartTime: 2019-10-10 11:20:30
	% EndTime: 2019-10-10 11:20:30
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
	% StartTime: 2019-10-10 11:20:30
	% EndTime: 2019-10-10 11:20:30
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (73->15), mult. (71->18), div. (0->0), fcn. (74->6), ass. (0->17)
	t32 = r_i_i_C(3) + qJ(4);
	t13 = qJ(2) + qJ(3);
	t10 = sin(t13);
	t11 = cos(t13);
	t27 = pkin(3) + r_i_i_C(1);
	t19 = t32 * t10 + t27 * t11;
	t31 = cos(qJ(2)) * pkin(2) + t19;
	t30 = t32 * t11;
	t28 = pkin(1) + t31;
	t15 = sin(qJ(1));
	t26 = t30 * t15;
	t16 = cos(qJ(1));
	t25 = t30 * t16;
	t23 = r_i_i_C(2) + pkin(8) + pkin(7);
	t20 = t27 * t10;
	t18 = -sin(qJ(2)) * pkin(2) - t20;
	t1 = [-t28 * t15 + t23 * t16, t18 * t16 + t25, -t16 * t20 + t25, t16 * t10, 0, 0; t23 * t15 + t28 * t16, t18 * t15 + t26, -t15 * t20 + t26, t15 * t10, 0, 0; 0, t31, t19, -t11, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:20:30
	% EndTime: 2019-10-10 11:20:30
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (92->19), mult. (87->20), div. (0->0), fcn. (92->6), ass. (0->19)
	t37 = r_i_i_C(1) + qJ(4);
	t15 = qJ(2) + qJ(3);
	t12 = sin(t15);
	t13 = cos(t15);
	t32 = -pkin(3) - pkin(4);
	t21 = t37 * t12 + (-r_i_i_C(2) - t32) * t13;
	t36 = cos(qJ(2)) * pkin(2) + t21;
	t35 = t37 * t13;
	t33 = pkin(1) + t36;
	t17 = sin(qJ(1));
	t30 = t17 * t12;
	t18 = cos(qJ(1));
	t29 = t18 * t12;
	t25 = r_i_i_C(2) * t30 + t35 * t17;
	t24 = r_i_i_C(2) * t29 + t35 * t18;
	t23 = -r_i_i_C(3) - qJ(5) + pkin(8) + pkin(7);
	t22 = t32 * t12;
	t20 = -sin(qJ(2)) * pkin(2) + t22;
	t1 = [-t33 * t17 + t23 * t18, t20 * t18 + t24, t18 * t22 + t24, t29, -t17, 0; t23 * t17 + t33 * t18, t20 * t17 + t25, t17 * t22 + t25, t30, t18, 0; 0, t36, t21, -t13, 0, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:20:30
	% EndTime: 2019-10-10 11:20:30
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (155->32), mult. (167->41), div. (0->0), fcn. (180->8), ass. (0->30)
	t26 = cos(qJ(6));
	t48 = pkin(5) + qJ(4);
	t50 = r_i_i_C(1) * t26 + t48;
	t22 = qJ(2) + qJ(3);
	t20 = cos(t22);
	t32 = pkin(3) + pkin(4) + pkin(9) + r_i_i_C(3);
	t49 = t32 * t20;
	t47 = t50 * t20;
	t19 = sin(t22);
	t21 = cos(qJ(2)) * pkin(2);
	t46 = t48 * t19 + pkin(1) + t21 + t49;
	t23 = sin(qJ(6));
	t43 = r_i_i_C(2) * t23;
	t25 = sin(qJ(1));
	t42 = t25 * t23;
	t41 = t25 * t26;
	t27 = cos(qJ(1));
	t40 = t27 * t23;
	t39 = t27 * t26;
	t36 = -qJ(5) + pkin(8) + pkin(7);
	t35 = t47 * t25;
	t34 = t47 * t27;
	t31 = -t32 * t19 - t20 * t43;
	t30 = t49 + (-t43 + t50) * t19;
	t29 = -sin(qJ(2)) * pkin(2) + t31;
	t4 = t19 * t39 - t42;
	t3 = -t19 * t40 - t41;
	t2 = -t19 * t41 - t40;
	t1 = t19 * t42 - t39;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t46 * t25 + t36 * t27, t29 * t27 + t34, t31 * t27 + t34, t27 * t19, -t25, t3 * r_i_i_C(1) - t4 * r_i_i_C(2); t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t36 * t25 + t46 * t27, t29 * t25 + t35, t31 * t25 + t35, t25 * t19, t27, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, t21 + t30, t30, -t20, 0, (r_i_i_C(1) * t23 + r_i_i_C(2) * t26) * t20;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end