% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPP7
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:51
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRPP7_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPP7_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP7_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_jacobia_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:51:06
	% EndTime: 2019-12-29 19:51:06
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:51:06
	% EndTime: 2019-12-29 19:51:06
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:51:05
	% EndTime: 2019-12-29 19:51:06
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
	% StartTime: 2019-12-29 19:51:11
	% EndTime: 2019-12-29 19:51:11
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (35->20), mult. (83->34), div. (0->0), fcn. (91->6), ass. (0->18)
	t16 = pkin(7) + r_i_i_C(3);
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
	t13 = [t10 * pkin(6) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t18 * t7, t11 * t10, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0; t7 * pkin(6) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t18 * t10, t11 * t7, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, t12 * t9 + t14, (-r_i_i_C(1) * t5 - r_i_i_C(2) * t8) * t6, 0, 0;];
	Ja_transl = t13;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:51:06
	% EndTime: 2019-12-29 19:51:06
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (61->22), mult. (142->33), div. (0->0), fcn. (163->6), ass. (0->22)
	t11 = cos(qJ(2));
	t20 = pkin(7) + r_i_i_C(2);
	t8 = sin(qJ(2));
	t15 = t20 * t8;
	t24 = t11 * pkin(2) + pkin(1) + t15;
	t10 = cos(qJ(3));
	t16 = r_i_i_C(3) + qJ(4);
	t21 = pkin(3) + r_i_i_C(1);
	t7 = sin(qJ(3));
	t23 = t21 * t10 + t16 * t7 + pkin(2);
	t9 = sin(qJ(1));
	t22 = t9 * t7;
	t12 = cos(qJ(1));
	t19 = t12 * t7;
	t18 = t9 * t10;
	t17 = t12 * t10;
	t13 = t20 * t11 - t23 * t8;
	t4 = t11 * t17 + t22;
	t3 = t11 * t19 - t18;
	t2 = t11 * t18 - t19;
	t1 = t11 * t22 + t17;
	t5 = [t12 * pkin(6) - t16 * t1 - t21 * t2 - t24 * t9, t13 * t12, t16 * t4 - t21 * t3, t3, 0; t9 * pkin(6) + t24 * t12 + t16 * t3 + t21 * t4, t13 * t9, -t21 * t1 + t16 * t2, t1, 0; 0, t23 * t11 + t15, (t16 * t10 - t21 * t7) * t8, t8 * t7, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:50:59
	% EndTime: 2019-12-29 19:50:59
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (80->24), mult. (179->36), div. (0->0), fcn. (207->6), ass. (0->21)
	t11 = cos(qJ(2));
	t16 = pkin(7) - r_i_i_C(3) - qJ(5);
	t8 = sin(qJ(2));
	t14 = t16 * t8;
	t23 = pkin(2) * t11 + pkin(1) + t14;
	t10 = cos(qJ(3));
	t17 = pkin(3) + pkin(4) + r_i_i_C(1);
	t18 = r_i_i_C(2) + qJ(4);
	t7 = sin(qJ(3));
	t22 = t17 * t10 + t18 * t7 + pkin(2);
	t9 = sin(qJ(1));
	t21 = t9 * t7;
	t20 = t9 * t10;
	t12 = cos(qJ(1));
	t19 = t11 * t12;
	t13 = t16 * t11 - t22 * t8;
	t4 = t10 * t19 + t21;
	t3 = t7 * t19 - t20;
	t2 = t11 * t20 - t12 * t7;
	t1 = t10 * t12 + t11 * t21;
	t5 = [pkin(6) * t12 - t18 * t1 - t17 * t2 - t23 * t9, t13 * t12, -t17 * t3 + t18 * t4, t3, -t12 * t8; t9 * pkin(6) + t23 * t12 + t17 * t4 + t18 * t3, t13 * t9, -t17 * t1 + t18 * t2, t1, -t9 * t8; 0, t22 * t11 + t14, (t18 * t10 - t17 * t7) * t8, t8 * t7, t11;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,5);
end