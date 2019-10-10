% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:49
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRP2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
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
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (8->5), mult. (6->6), div. (0->0), fcn. (6->4), ass. (0->4)
	t3 = qJ(1) + pkin(9);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [-t1 * r_i_i_C(1) - t2 * r_i_i_C(2) - sin(qJ(1)) * pkin(1), 0, 0, 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + cos(qJ(1)) * pkin(1), 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (23->9), mult. (16->8), div. (0->0), fcn. (18->6), ass. (0->6)
	t7 = r_i_i_C(3) + qJ(3);
	t6 = r_i_i_C(1) * cos(pkin(10)) - r_i_i_C(2) * sin(pkin(10)) + pkin(2);
	t3 = qJ(1) + pkin(9);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [-sin(qJ(1)) * pkin(1) + t7 * t2 - t6 * t1, 0, t1, 0, 0, 0; cos(qJ(1)) * pkin(1) + t7 * t1 + t6 * t2, 0, -t2, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (44->12), mult. (28->13), div. (0->0), fcn. (30->7), ass. (0->11)
	t12 = r_i_i_C(3) + pkin(7) + qJ(3);
	t6 = pkin(10) + qJ(4);
	t2 = sin(t6);
	t4 = cos(t6);
	t11 = r_i_i_C(1) * t4 - r_i_i_C(2) * t2;
	t10 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t4;
	t9 = cos(pkin(10)) * pkin(3) + pkin(2) + t11;
	t7 = qJ(1) + pkin(9);
	t5 = cos(t7);
	t3 = sin(t7);
	t1 = [-sin(qJ(1)) * pkin(1) + t12 * t5 - t9 * t3, 0, t3, t10 * t5, 0, 0; cos(qJ(1)) * pkin(1) + t12 * t3 + t9 * t5, 0, -t5, t10 * t3, 0, 0; 0, 1, 0, t11, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (111->27), mult. (87->35), div. (0->0), fcn. (97->9), ass. (0->23)
	t23 = pkin(8) + r_i_i_C(3);
	t10 = pkin(10) + qJ(4);
	t6 = sin(t10);
	t18 = t23 * t6;
	t8 = cos(t10);
	t24 = t18 + t8 * pkin(4) + cos(pkin(10)) * pkin(3) + pkin(2);
	t14 = cos(qJ(5));
	t11 = qJ(1) + pkin(9);
	t7 = sin(t11);
	t22 = t14 * t7;
	t9 = cos(t11);
	t21 = t14 * t9;
	t13 = sin(qJ(5));
	t20 = t7 * t13;
	t19 = t9 * t13;
	t16 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + pkin(4);
	t15 = -t16 * t6 + t23 * t8;
	t12 = -pkin(7) - qJ(3);
	t4 = t8 * t21 + t20;
	t3 = -t8 * t19 + t22;
	t2 = -t8 * t22 + t19;
	t1 = t8 * t20 + t21;
	t5 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t9 * t12 - t24 * t7, 0, t7, t15 * t9, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t7 * t12 + t24 * t9, 0, -t9, t15 * t7, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, 1, 0, t16 * t8 + t18, (-r_i_i_C(1) * t13 - r_i_i_C(2) * t14) * t6, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (180->29), mult. (146->36), div. (0->0), fcn. (169->9), ass. (0->25)
	t12 = pkin(10) + qJ(4);
	t10 = cos(t12);
	t25 = pkin(8) + r_i_i_C(2);
	t8 = sin(t12);
	t19 = t25 * t8;
	t28 = t19 + t10 * pkin(4) + cos(pkin(10)) * pkin(3) + pkin(2);
	t15 = sin(qJ(5));
	t16 = cos(qJ(5));
	t20 = r_i_i_C(3) + qJ(6);
	t26 = pkin(5) + r_i_i_C(1);
	t27 = t20 * t15 + t26 * t16 + pkin(4);
	t13 = qJ(1) + pkin(9);
	t9 = sin(t13);
	t24 = t9 * t15;
	t23 = t9 * t16;
	t11 = cos(t13);
	t22 = t11 * t15;
	t21 = t11 * t16;
	t17 = t25 * t10 - t27 * t8;
	t14 = -pkin(7) - qJ(3);
	t4 = t10 * t21 + t24;
	t3 = t10 * t22 - t23;
	t2 = t10 * t23 - t22;
	t1 = t10 * t24 + t21;
	t5 = [-sin(qJ(1)) * pkin(1) - t11 * t14 - t26 * t2 - t20 * t1 - t28 * t9, 0, t9, t17 * t11, t20 * t4 - t26 * t3, t3; cos(qJ(1)) * pkin(1) - t9 * t14 + t26 * t4 + t20 * t3 + t28 * t11, 0, -t11, t17 * t9, -t26 * t1 + t20 * t2, t1; 0, 1, 0, t27 * t10 + t19, (-t26 * t15 + t20 * t16) * t8, t8 * t15;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end