% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:19
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPPR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
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
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.11s
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
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (53->13), mult. (51->14), div. (0->0), fcn. (56->6), ass. (0->12)
	t12 = r_i_i_C(3) + qJ(4);
	t14 = pkin(3) - r_i_i_C(2);
	t5 = qJ(2) + pkin(9);
	t2 = sin(t5);
	t3 = cos(t5);
	t16 = cos(qJ(2)) * pkin(2) + t12 * t2 + t14 * t3;
	t15 = pkin(1) + t16;
	t13 = r_i_i_C(1) + qJ(3) + pkin(7);
	t10 = -sin(qJ(2)) * pkin(2) + t12 * t3 - t14 * t2;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [t13 * t9 - t15 * t8, t10 * t9, t8, t9 * t2, 0, 0; t13 * t8 + t15 * t9, t10 * t8, -t9, t8 * t2, 0, 0; 0, t16, 0, -t3, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (80->17), mult. (85->20), div. (0->0), fcn. (97->8), ass. (0->14)
	t6 = sin(pkin(10));
	t7 = cos(pkin(10));
	t15 = r_i_i_C(1) * t6 + r_i_i_C(2) * t7 + qJ(4);
	t16 = pkin(3) + r_i_i_C(3) + qJ(5);
	t5 = qJ(2) + pkin(9);
	t2 = sin(t5);
	t3 = cos(t5);
	t18 = cos(qJ(2)) * pkin(2) + t15 * t2 + t16 * t3;
	t17 = pkin(1) + t18;
	t14 = r_i_i_C(1) * t7 - t6 * r_i_i_C(2) + pkin(4) + pkin(7) + qJ(3);
	t12 = -sin(qJ(2)) * pkin(2) + t15 * t3 - t16 * t2;
	t11 = cos(qJ(1));
	t10 = sin(qJ(1));
	t1 = [-t17 * t10 + t14 * t11, t12 * t11, t10, t11 * t2, t11 * t3, 0; t14 * t10 + t17 * t11, t12 * t10, -t11, t10 * t2, t10 * t3, 0; 0, t18, 0, -t3, t2, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (134->30), mult. (121->40), div. (0->0), fcn. (137->10), ass. (0->24)
	t13 = qJ(2) + pkin(9);
	t10 = cos(t13);
	t24 = pkin(3) + r_i_i_C(3) + pkin(8) + qJ(5);
	t32 = cos(qJ(2)) * pkin(2) + t24 * t10;
	t23 = pkin(5) * sin(pkin(10)) + qJ(4);
	t8 = sin(t13);
	t31 = -t23 * t8 - pkin(1) - t32;
	t18 = sin(qJ(1));
	t29 = t18 * t8;
	t12 = pkin(10) + qJ(6);
	t9 = cos(t12);
	t28 = t18 * t9;
	t19 = cos(qJ(1));
	t27 = t19 * t8;
	t26 = t19 * t9;
	t25 = qJ(3) + pkin(7) + cos(pkin(10)) * pkin(5) + pkin(4);
	t7 = sin(t12);
	t21 = r_i_i_C(1) * t7 + r_i_i_C(2) * t9 + t23;
	t20 = -sin(qJ(2)) * pkin(2) + t21 * t10 - t24 * t8;
	t4 = -t7 * t29 + t26;
	t3 = t19 * t7 + t8 * t28;
	t2 = t7 * t27 + t28;
	t1 = -t18 * t7 + t8 * t26;
	t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t31 * t18 + t25 * t19, t20 * t19, t18, t27, t19 * t10, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t25 * t18 - t31 * t19, t20 * t18, -t19, t29, t18 * t10, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, t21 * t8 + t32, 0, -t10, t8, (-r_i_i_C(1) * t9 + r_i_i_C(2) * t7) * t10;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end