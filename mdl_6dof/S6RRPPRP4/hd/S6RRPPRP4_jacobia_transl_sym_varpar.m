% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:32
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRP4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRP4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:32:13
	% EndTime: 2019-10-10 09:32:13
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:32:13
	% EndTime: 2019-10-10 09:32:13
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
	% StartTime: 2019-10-10 09:32:13
	% EndTime: 2019-10-10 09:32:13
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
	% StartTime: 2019-10-10 09:32:13
	% EndTime: 2019-10-10 09:32:13
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (29->11), mult. (65->16), div. (0->0), fcn. (72->6), ass. (0->13)
	t1 = sin(pkin(9));
	t2 = cos(pkin(9));
	t10 = r_i_i_C(1) * t2 - r_i_i_C(2) * t1 + pkin(2);
	t11 = r_i_i_C(3) + qJ(3);
	t3 = sin(qJ(2));
	t5 = cos(qJ(2));
	t8 = t10 * t5 + t11 * t3;
	t12 = pkin(1) + t8;
	t9 = t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + pkin(7);
	t7 = -t10 * t3 + t11 * t5;
	t6 = cos(qJ(1));
	t4 = sin(qJ(1));
	t13 = [-t12 * t4 + t9 * t6, t7 * t6, t6 * t3, 0, 0, 0; t12 * t6 + t9 * t4, t7 * t4, t4 * t3, 0, 0, 0; 0, t8, -t5, 0, 0, 0;];
	Ja_transl = t13;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:32:13
	% EndTime: 2019-10-10 09:32:13
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (45->20), mult. (104->30), div. (0->0), fcn. (120->6), ass. (0->18)
	t15 = r_i_i_C(2) + qJ(3);
	t7 = sin(qJ(2));
	t12 = t15 * t7;
	t9 = cos(qJ(2));
	t20 = t9 * pkin(2) + pkin(1) + t12;
	t14 = r_i_i_C(3) + qJ(4);
	t17 = pkin(3) + r_i_i_C(1);
	t5 = sin(pkin(9));
	t6 = cos(pkin(9));
	t19 = t14 * t5 + t17 * t6 + pkin(2);
	t8 = sin(qJ(1));
	t18 = t8 * t9;
	t10 = cos(qJ(1));
	t16 = t10 * t9;
	t11 = t15 * t9 - t19 * t7;
	t3 = t5 * t16 - t8 * t6;
	t1 = t10 * t6 + t5 * t18;
	t2 = [pkin(7) * t10 + t17 * (t10 * t5 - t6 * t18) - t14 * t1 - t20 * t8, t11 * t10, t10 * t7, t3, 0, 0; t8 * pkin(7) + t17 * (t6 * t16 + t8 * t5) + t14 * t3 + t20 * t10, t11 * t8, t8 * t7, t1, 0, 0; 0, t19 * t9 + t12, -t9, t7 * t5, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:32:13
	% EndTime: 2019-10-10 09:32:13
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (88->34), mult. (217->55), div. (0->0), fcn. (263->8), ass. (0->25)
	t14 = cos(qJ(2));
	t11 = sin(qJ(2));
	t21 = -r_i_i_C(3) - pkin(8) + qJ(3);
	t19 = t21 * t11;
	t26 = t14 * pkin(2) + pkin(1) + t19;
	t10 = sin(qJ(5));
	t13 = cos(qJ(5));
	t24 = pkin(3) + pkin(4);
	t17 = t13 * r_i_i_C(1) - t10 * r_i_i_C(2) + t24;
	t18 = t10 * r_i_i_C(1) + t13 * r_i_i_C(2) + qJ(4);
	t8 = sin(pkin(9));
	t9 = cos(pkin(9));
	t25 = t17 * t9 + t18 * t8 + pkin(2);
	t12 = sin(qJ(1));
	t23 = t12 * t14;
	t15 = cos(qJ(1));
	t22 = t14 * t15;
	t16 = -t25 * t11 + t21 * t14;
	t6 = t12 * t8 + t9 * t22;
	t5 = -t12 * t9 + t8 * t22;
	t4 = -t15 * t8 + t9 * t23;
	t3 = t15 * t9 + t8 * t23;
	t2 = t5 * t10 + t13 * t6;
	t1 = -t6 * t10 + t13 * t5;
	t7 = [t15 * pkin(7) - t26 * t12 - t17 * t4 - t18 * t3, t16 * t15, t15 * t11, t5, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t12 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t5 * qJ(4) + t26 * t15 + t24 * t6, t16 * t12, t12 * t11, t3, (-t4 * t10 + t13 * t3) * r_i_i_C(1) + (-t3 * t10 - t13 * t4) * r_i_i_C(2), 0; 0, t25 * t14 + t19, -t14, t11 * t8, ((-t10 * t9 + t13 * t8) * r_i_i_C(1) + (-t10 * t8 - t13 * t9) * r_i_i_C(2)) * t11, 0;];
	Ja_transl = t7;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:32:13
	% EndTime: 2019-10-10 09:32:13
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (143->38), mult. (357->56), div. (0->0), fcn. (445->8), ass. (0->32)
	t23 = sin(pkin(9));
	t24 = cos(pkin(9));
	t28 = cos(qJ(5));
	t25 = sin(qJ(5));
	t34 = t23 * t25 + t24 * t28;
	t42 = r_i_i_C(3) + qJ(6);
	t47 = t24 * t25;
	t48 = pkin(5) + r_i_i_C(1);
	t49 = pkin(3) + pkin(4);
	t55 = -qJ(4) * t23 - t49 * t24 - t48 * t34 - pkin(2) - t42 * (-t23 * t28 + t47);
	t29 = cos(qJ(2));
	t26 = sin(qJ(2));
	t41 = pkin(8) + r_i_i_C(2) - qJ(3);
	t38 = t41 * t26;
	t54 = -t29 * pkin(2) - pkin(1) + t38;
	t51 = t55 * t26 - t41 * t29;
	t46 = t26 * t23;
	t27 = sin(qJ(1));
	t45 = t27 * t29;
	t30 = cos(qJ(1));
	t44 = t30 * t23;
	t43 = t30 * t24;
	t17 = t23 * t45 + t43;
	t18 = t24 * t45 - t44;
	t37 = t17 * t28 - t18 * t25;
	t36 = t17 * t25 + t18 * t28;
	t20 = t27 * t23 + t29 * t43;
	t19 = -t27 * t24 + t29 * t44;
	t13 = t26 * t47 - t28 * t46;
	t6 = t19 * t25 + t20 * t28;
	t5 = -t19 * t28 + t20 * t25;
	t1 = [t30 * pkin(7) - t17 * qJ(4) - t49 * t18 + t54 * t27 - t48 * t36 + t42 * t37, t51 * t30, t30 * t26, t19, t42 * t6 - t48 * t5, t5; t27 * pkin(7) + t19 * qJ(4) + t49 * t20 - t54 * t30 + t42 * t5 + t48 * t6, t51 * t27, t27 * t26, t17, t42 * t36 + t37 * t48, -t37; 0, -t55 * t29 - t38, -t29, t46, t42 * t34 * t26 - t48 * t13, t13;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end