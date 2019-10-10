% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:57
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPP2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPP2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:57:25
	% EndTime: 2019-10-10 09:57:25
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:57:25
	% EndTime: 2019-10-10 09:57:25
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
	% StartTime: 2019-10-10 09:57:25
	% EndTime: 2019-10-10 09:57:25
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
	% StartTime: 2019-10-10 09:57:25
	% EndTime: 2019-10-10 09:57:25
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (29->10), mult. (31->12), div. (0->0), fcn. (33->6), ass. (0->10)
	t5 = qJ(2) + pkin(9);
	t2 = sin(t5);
	t3 = cos(t5);
	t14 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2) + cos(qJ(2)) * pkin(2);
	t13 = r_i_i_C(3) + qJ(3) + pkin(7);
	t11 = pkin(1) + t14;
	t10 = -sin(qJ(2)) * pkin(2) - r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t11 * t8 + t13 * t9, t10 * t9, t8, 0, 0, 0; t11 * t9 + t13 * t8, t10 * t8, -t9, 0, 0, 0; 0, t14, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:57:25
	% EndTime: 2019-10-10 09:57:25
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (74->25), mult. (90->34), div. (0->0), fcn. (100->8), ass. (0->22)
	t24 = pkin(8) + r_i_i_C(3);
	t9 = qJ(2) + pkin(9);
	t6 = sin(t9);
	t26 = t24 * t6 + cos(qJ(2)) * pkin(2);
	t7 = cos(t9);
	t25 = t7 * pkin(3) + pkin(1) + t26;
	t11 = sin(qJ(4));
	t15 = cos(qJ(1));
	t23 = t11 * t15;
	t13 = sin(qJ(1));
	t22 = t13 * t11;
	t14 = cos(qJ(4));
	t21 = t13 * t14;
	t20 = t14 * t15;
	t17 = r_i_i_C(1) * t14 - r_i_i_C(2) * t11 + pkin(3);
	t16 = -sin(qJ(2)) * pkin(2) - t17 * t6 + t24 * t7;
	t10 = -qJ(3) - pkin(7);
	t4 = t7 * t20 + t22;
	t3 = -t7 * t23 + t21;
	t2 = -t7 * t21 + t23;
	t1 = t7 * t22 + t20;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t15 * t10 - t25 * t13, t16 * t15, t13, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t13 * t10 + t25 * t15, t16 * t13, -t15, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, t17 * t7 + t26, 0, (-r_i_i_C(1) * t11 - r_i_i_C(2) * t14) * t6, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:57:25
	% EndTime: 2019-10-10 09:57:25
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (119->27), mult. (149->35), div. (0->0), fcn. (172->8), ass. (0->24)
	t26 = pkin(8) + r_i_i_C(2);
	t11 = qJ(2) + pkin(9);
	t8 = sin(t11);
	t30 = cos(qJ(2)) * pkin(2) + t26 * t8;
	t9 = cos(t11);
	t29 = pkin(3) * t9 + pkin(1) + t30;
	t13 = sin(qJ(4));
	t16 = cos(qJ(4));
	t21 = r_i_i_C(3) + qJ(5);
	t27 = pkin(4) + r_i_i_C(1);
	t28 = t21 * t13 + t27 * t16 + pkin(3);
	t15 = sin(qJ(1));
	t25 = t15 * t13;
	t24 = t15 * t16;
	t17 = cos(qJ(1));
	t23 = t16 * t17;
	t22 = t17 * t13;
	t18 = -sin(qJ(2)) * pkin(2) + t26 * t9 - t28 * t8;
	t12 = -qJ(3) - pkin(7);
	t4 = t9 * t23 + t25;
	t3 = t9 * t22 - t24;
	t2 = t9 * t24 - t22;
	t1 = t9 * t25 + t23;
	t5 = [-t21 * t1 - t12 * t17 - t29 * t15 - t27 * t2, t18 * t17, t15, t21 * t4 - t27 * t3, t3, 0; -t15 * t12 + t29 * t17 + t21 * t3 + t27 * t4, t18 * t15, -t17, -t27 * t1 + t21 * t2, t1, 0; 0, t28 * t9 + t30, 0, (-t27 * t13 + t21 * t16) * t8, t8 * t13, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:57:25
	% EndTime: 2019-10-10 09:57:25
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (154->29), mult. (186->37), div. (0->0), fcn. (216->8), ass. (0->24)
	t21 = pkin(8) - r_i_i_C(3) - qJ(6);
	t11 = qJ(2) + pkin(9);
	t8 = sin(t11);
	t30 = cos(qJ(2)) * pkin(2) + t21 * t8;
	t9 = cos(t11);
	t29 = t9 * pkin(3) + pkin(1) + t30;
	t13 = sin(qJ(4));
	t16 = cos(qJ(4));
	t22 = pkin(4) + pkin(5) + r_i_i_C(1);
	t23 = r_i_i_C(2) + qJ(5);
	t28 = t23 * t13 + t22 * t16 + pkin(3);
	t15 = sin(qJ(1));
	t27 = t15 * t13;
	t26 = t15 * t16;
	t17 = cos(qJ(1));
	t25 = t17 * t13;
	t24 = t17 * t16;
	t18 = -sin(qJ(2)) * pkin(2) + t21 * t9 - t28 * t8;
	t12 = -qJ(3) - pkin(7);
	t4 = t9 * t24 + t27;
	t3 = t9 * t25 - t26;
	t2 = t9 * t26 - t25;
	t1 = t9 * t27 + t24;
	t5 = [-t23 * t1 - t17 * t12 - t29 * t15 - t22 * t2, t18 * t17, t15, -t22 * t3 + t23 * t4, t3, -t17 * t8; -t15 * t12 + t29 * t17 + t22 * t4 + t23 * t3, t18 * t15, -t17, -t22 * t1 + t23 * t2, t1, -t15 * t8; 0, t28 * t9 + t30, 0, (-t22 * t13 + t23 * t16) * t8, t8 * t13, t9;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end