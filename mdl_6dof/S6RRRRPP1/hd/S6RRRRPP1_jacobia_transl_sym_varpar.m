% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:20
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPP1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
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
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
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
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
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
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (100->28), mult. (121->39), div. (0->0), fcn. (129->8), ass. (0->29)
	t19 = qJ(2) + qJ(3);
	t16 = sin(t19);
	t17 = cos(t19);
	t20 = sin(qJ(4));
	t39 = r_i_i_C(2) * t20;
	t45 = pkin(9) + r_i_i_C(3);
	t46 = t16 * t39 + t17 * t45;
	t43 = t17 * pkin(3) + t45 * t16;
	t18 = cos(qJ(2)) * pkin(2);
	t42 = pkin(1) + t18 + t43;
	t23 = cos(qJ(4));
	t40 = r_i_i_C(1) * t23;
	t22 = sin(qJ(1));
	t36 = t20 * t22;
	t24 = cos(qJ(1));
	t35 = t20 * t24;
	t34 = t22 * t23;
	t33 = t23 * t24;
	t32 = t46 * t22;
	t30 = t46 * t24;
	t28 = (-pkin(3) - t40) * t16;
	t27 = (-t39 + t40) * t17 + t43;
	t26 = -sin(qJ(2)) * pkin(2) + t28;
	t25 = -pkin(8) - pkin(7);
	t4 = t17 * t33 + t36;
	t3 = -t17 * t35 + t34;
	t2 = -t17 * t34 + t35;
	t1 = t17 * t36 + t33;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t42 * t22 - t24 * t25, t26 * t24 + t30, t24 * t28 + t30, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t22 * t25 + t42 * t24, t26 * t22 + t32, t22 * t28 + t32, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, t18 + t27, t27, (-r_i_i_C(1) * t20 - r_i_i_C(2) * t23) * t16, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (160->39), mult. (147->53), div. (0->0), fcn. (158->10), ass. (0->36)
	t25 = cos(qJ(4));
	t12 = t25 * pkin(4) + pkin(3);
	t20 = qJ(2) + qJ(3);
	t16 = sin(t20);
	t17 = cos(t20);
	t21 = -qJ(5) - pkin(9);
	t47 = t17 * t12 + (r_i_i_C(3) - t21) * t16;
	t18 = cos(qJ(2)) * pkin(2);
	t46 = pkin(1) + t18 + t47;
	t24 = sin(qJ(1));
	t19 = qJ(4) + pkin(10);
	t14 = sin(t19);
	t42 = r_i_i_C(2) * t14;
	t33 = t16 * t42;
	t39 = t17 * t24;
	t45 = r_i_i_C(3) * t39 + t24 * t33;
	t22 = sin(qJ(4));
	t44 = pkin(4) * t22;
	t15 = cos(t19);
	t43 = r_i_i_C(1) * t15;
	t26 = cos(qJ(1));
	t38 = t17 * t26;
	t40 = r_i_i_C(3) * t38 + t26 * t33;
	t37 = t24 * t14;
	t36 = t24 * t15;
	t35 = t26 * t14;
	t34 = t26 * t15;
	t32 = pkin(8) + pkin(7) + t44;
	t30 = -t17 * t21 + (-t12 - t43) * t16;
	t29 = (-t42 + t43) * t17 + t47;
	t28 = -sin(qJ(2)) * pkin(2) + t30;
	t4 = t17 * t34 + t37;
	t3 = -t17 * t35 + t36;
	t2 = -t17 * t36 + t35;
	t1 = t17 * t37 + t34;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t46 * t24 + t32 * t26, t28 * t26 + t40, t30 * t26 + t40, t3 * r_i_i_C(1) - t4 * r_i_i_C(2) + (-t22 * t38 + t24 * t25) * pkin(4), t26 * t16, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t32 * t24 + t46 * t26, t28 * t24 + t45, t30 * t24 + t45, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + (-t22 * t39 - t25 * t26) * pkin(4), t24 * t16, 0; 0, t18 + t29, t29, (-r_i_i_C(1) * t14 - r_i_i_C(2) * t15 - t44) * t16, -t17, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (252->39), mult. (222->51), div. (0->0), fcn. (246->10), ass. (0->36)
	t22 = qJ(4) + pkin(10);
	t17 = sin(t22);
	t18 = cos(t22);
	t36 = r_i_i_C(3) + qJ(6);
	t47 = pkin(5) + r_i_i_C(1);
	t51 = t36 * t17 + t47 * t18;
	t28 = cos(qJ(4));
	t15 = t28 * pkin(4) + pkin(3);
	t23 = qJ(2) + qJ(3);
	t19 = sin(t23);
	t20 = cos(t23);
	t24 = -qJ(5) - pkin(9);
	t49 = t20 * t15 + (r_i_i_C(2) - t24) * t19;
	t21 = cos(qJ(2)) * pkin(2);
	t48 = pkin(1) + t21 + t49;
	t25 = sin(qJ(4));
	t46 = pkin(4) * t25;
	t27 = sin(qJ(1));
	t42 = t20 * t27;
	t29 = cos(qJ(1));
	t41 = t20 * t29;
	t40 = t27 * t17;
	t39 = t27 * t18;
	t38 = t29 * t17;
	t37 = t29 * t18;
	t35 = pkin(8) + pkin(7) + t46;
	t33 = t51 * t20 + t49;
	t32 = -t20 * t24 + (-t15 - t51) * t19;
	t31 = -sin(qJ(2)) * pkin(2) + t32;
	t11 = r_i_i_C(2) * t41;
	t10 = r_i_i_C(2) * t42;
	t4 = t20 * t37 + t40;
	t3 = t20 * t38 - t39;
	t2 = t20 * t39 - t38;
	t1 = t20 * t40 + t37;
	t5 = [-t36 * t1 - t47 * t2 - t48 * t27 + t35 * t29, t31 * t29 + t11, t32 * t29 + t11, t36 * t4 - t47 * t3 + (-t25 * t41 + t27 * t28) * pkin(4), t29 * t19, t3; t35 * t27 + t48 * t29 + t36 * t3 + t47 * t4, t31 * t27 + t10, t32 * t27 + t10, t36 * t2 - t47 * t1 + (-t25 * t42 - t28 * t29) * pkin(4), t27 * t19, t1; 0, t21 + t33, t33, (-t47 * t17 + t36 * t18 - t46) * t19, -t20, t19 * t17;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end