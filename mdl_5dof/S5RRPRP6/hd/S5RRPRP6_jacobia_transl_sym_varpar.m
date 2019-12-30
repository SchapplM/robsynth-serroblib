% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRP6
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:47
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRPRP6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRP6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:47:50
	% EndTime: 2019-12-29 18:47:50
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:47:50
	% EndTime: 2019-12-29 18:47:50
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:47:50
	% EndTime: 2019-12-29 18:47:50
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (11->6), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->9)
	t8 = pkin(6) + r_i_i_C(3);
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = t3 * r_i_i_C(1) - t1 * r_i_i_C(2);
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t5 * t2 + t8 * t4, t6 * t4, 0, 0, 0; t8 * t2 + t5 * t4, t6 * t2, 0, 0, 0; 0, t7, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:47:50
	% EndTime: 2019-12-29 18:47:50
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (29->10), mult. (31->12), div. (0->0), fcn. (33->6), ass. (0->10)
	t5 = qJ(2) + pkin(8);
	t2 = sin(t5);
	t3 = cos(t5);
	t14 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2) + cos(qJ(2)) * pkin(2);
	t13 = r_i_i_C(3) + qJ(3) + pkin(6);
	t11 = pkin(1) + t14;
	t10 = -sin(qJ(2)) * pkin(2) - r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t11 * t8 + t13 * t9, t10 * t9, t8, 0, 0; t11 * t9 + t13 * t8, t10 * t8, -t9, 0, 0; 0, t14, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:47:50
	% EndTime: 2019-12-29 18:47:50
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (74->25), mult. (90->34), div. (0->0), fcn. (100->8), ass. (0->22)
	t24 = pkin(7) + r_i_i_C(3);
	t9 = qJ(2) + pkin(8);
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
	t10 = -qJ(3) - pkin(6);
	t4 = t7 * t20 + t22;
	t3 = -t7 * t23 + t21;
	t2 = -t7 * t21 + t23;
	t1 = t7 * t22 + t20;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t10 * t15 - t25 * t13, t16 * t15, t13, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t13 * t10 + t25 * t15, t16 * t13, -t15, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, t17 * t7 + t26, 0, (-r_i_i_C(1) * t11 - r_i_i_C(2) * t14) * t6, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:47:43
	% EndTime: 2019-12-29 18:47:43
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (98->29), mult. (113->38), div. (0->0), fcn. (126->8), ass. (0->24)
	t27 = r_i_i_C(3) + qJ(5) + pkin(7);
	t10 = qJ(2) + pkin(8);
	t7 = sin(t10);
	t30 = t27 * t7 + cos(qJ(2)) * pkin(2);
	t16 = cos(qJ(4));
	t5 = pkin(4) * t16 + pkin(3);
	t8 = cos(t10);
	t29 = t5 * t8 + pkin(1) + t30;
	t28 = pkin(4) + r_i_i_C(1);
	t13 = sin(qJ(4));
	t17 = cos(qJ(1));
	t26 = t13 * t17;
	t15 = sin(qJ(1));
	t25 = t15 * t13;
	t24 = t15 * t16;
	t23 = t16 * t17;
	t20 = pkin(4) * t13 + pkin(6) + qJ(3);
	t19 = r_i_i_C(1) * t16 - r_i_i_C(2) * t13 + t5;
	t3 = -t8 * t26 + t24;
	t1 = t8 * t25 + t23;
	t18 = -sin(qJ(2)) * pkin(2) - t19 * t7 + t27 * t8;
	t4 = t8 * t23 + t25;
	t2 = -t8 * t24 + t26;
	t6 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t29 * t15 + t20 * t17, t18 * t17, t15, -t4 * r_i_i_C(2) + t28 * t3, t17 * t7; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t20 * t15 + t29 * t17, t18 * t15, -t17, t2 * r_i_i_C(2) - t28 * t1, t15 * t7; 0, t19 * t8 + t30, 0, (-r_i_i_C(2) * t16 - t28 * t13) * t7, -t8;];
	Ja_transl = t6;
else
	Ja_transl=NaN(3,5);
end