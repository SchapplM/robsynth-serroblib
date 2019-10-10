% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRPR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:39:24
	% EndTime: 2019-10-09 23:39:24
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:39:24
	% EndTime: 2019-10-09 23:39:24
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
	% StartTime: 2019-10-09 23:39:24
	% EndTime: 2019-10-09 23:39:24
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (7->5), mult. (8->4), div. (0->0), fcn. (10->2), ass. (0->5)
	t4 = pkin(1) + r_i_i_C(1);
	t3 = r_i_i_C(3) + qJ(2);
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t5 = [-t1 * t4 + t2 * t3, t1, 0, 0, 0, 0; t1 * t3 + t2 * t4, -t2, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:39:24
	% EndTime: 2019-10-09 23:39:24
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (13->10), mult. (18->12), div. (0->0), fcn. (24->4), ass. (0->8)
	t7 = pkin(1) + pkin(2);
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t4 = cos(pkin(9));
	t3 = sin(pkin(9));
	t2 = t3 * t5 + t4 * t6;
	t1 = t3 * t6 - t4 * t5;
	t8 = [t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + t6 * qJ(2) - t5 * t7, t5, 0, 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + t5 * qJ(2) + t6 * t7, -t6, 0, 0, 0, 0; 0, 0, -1, 0, 0, 0;];
	Ja_transl = t8;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:39:24
	% EndTime: 2019-10-09 23:39:24
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (32->14), mult. (58->18), div. (0->0), fcn. (74->6), ass. (0->14)
	t16 = pkin(1) + pkin(2);
	t15 = pkin(7) + r_i_i_C(3);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t12 = cos(pkin(9));
	t11 = sin(pkin(9));
	t6 = sin(qJ(4));
	t7 = cos(qJ(4));
	t10 = -t7 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t9 = r_i_i_C(1) * t6 + r_i_i_C(2) * t7;
	t8 = pkin(3) - t10;
	t2 = t14 * t11 - t13 * t12;
	t1 = -t13 * t11 - t14 * t12;
	t3 = [t14 * qJ(2) + t15 * t1 - t16 * t13 + t8 * t2, t13, 0, t9 * t1, 0, 0; t13 * qJ(2) - t8 * t1 + t16 * t14 + t15 * t2, -t14, 0, t9 * t2, 0, 0; 0, 0, -1, t10, 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:39:24
	% EndTime: 2019-10-09 23:39:24
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (53->18), mult. (73->20), div. (0->0), fcn. (95->8), ass. (0->15)
	t23 = pkin(1) + pkin(2);
	t9 = qJ(4) + pkin(10);
	t7 = sin(t9);
	t8 = cos(t9);
	t22 = -t8 * r_i_i_C(1) + t7 * r_i_i_C(2) - cos(qJ(4)) * pkin(4);
	t20 = r_i_i_C(3) + qJ(5) + pkin(7);
	t19 = cos(qJ(1));
	t18 = sin(qJ(1));
	t17 = cos(pkin(9));
	t16 = sin(pkin(9));
	t14 = pkin(3) - t22;
	t13 = sin(qJ(4)) * pkin(4) + r_i_i_C(1) * t7 + r_i_i_C(2) * t8;
	t2 = t19 * t16 - t18 * t17;
	t1 = -t18 * t16 - t19 * t17;
	t3 = [t19 * qJ(2) + t20 * t1 + t14 * t2 - t23 * t18, t18, 0, t13 * t1, t2, 0; t18 * qJ(2) - t14 * t1 + t23 * t19 + t20 * t2, -t19, 0, t13 * t2, -t1, 0; 0, 0, -1, t22, 0, 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:39:24
	% EndTime: 2019-10-09 23:39:24
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (120->34), mult. (176->42), div. (0->0), fcn. (228->10), ass. (0->26)
	t11 = qJ(4) + pkin(10);
	t10 = cos(t11);
	t13 = sin(qJ(6));
	t15 = cos(qJ(6));
	t18 = t15 * r_i_i_C(1) - t13 * r_i_i_C(2) + pkin(5);
	t28 = pkin(8) + r_i_i_C(3);
	t9 = sin(t11);
	t20 = t28 * t9;
	t31 = -t18 * t10 - t20;
	t30 = pkin(1) + pkin(2);
	t27 = cos(qJ(4)) * pkin(4);
	t26 = cos(qJ(1));
	t25 = sin(qJ(1));
	t24 = t10 * t13;
	t23 = t10 * t15;
	t22 = cos(pkin(9));
	t21 = sin(pkin(9));
	t19 = t13 * r_i_i_C(1) + t15 * r_i_i_C(2);
	t17 = sin(qJ(4)) * pkin(4) - t28 * t10 + t18 * t9;
	t12 = -qJ(5) - pkin(7);
	t8 = pkin(3) + t27;
	t4 = t26 * t21 - t25 * t22;
	t3 = -t25 * t21 - t26 * t22;
	t2 = t4 * t13 - t3 * t23;
	t1 = t4 * t15 + t3 * t24;
	t5 = [t26 * qJ(2) + (-t12 + t19) * t3 + (t8 - t31) * t4 - t30 * t25, t25, 0, t17 * t3, t4, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t25 * qJ(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t4 * t12 + (-t10 * pkin(5) - t20 - t8) * t3 + t30 * t26, -t26, 0, t17 * t4, -t3, (-t3 * t15 + t4 * t24) * r_i_i_C(1) + (t3 * t13 + t4 * t23) * r_i_i_C(2); 0, 0, -1, -t27 + t31, 0, t19 * t9;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end