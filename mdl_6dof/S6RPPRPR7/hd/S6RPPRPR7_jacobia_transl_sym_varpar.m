% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRPR7
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
% Datum: 2019-10-09 23:44
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRPR7_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR7_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR7_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
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
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (7->5), mult. (8->4), div. (0->0), fcn. (10->2), ass. (0->5)
	t4 = pkin(1) - r_i_i_C(2);
	t3 = r_i_i_C(3) + qJ(2);
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t5 = [-t1 * t4 + t2 * t3, t1, 0, 0, 0, 0; t1 * t3 + t2 * t4, -t2, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (11->6), mult. (16->6), div. (0->0), fcn. (20->4), ass. (0->5)
	t6 = pkin(1) + r_i_i_C(3) + qJ(3);
	t5 = r_i_i_C(1) * sin(pkin(9)) + r_i_i_C(2) * cos(pkin(9)) + qJ(2);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t6 * t3 + t5 * t4, t3, t4, 0, 0, 0; t5 * t3 + t6 * t4, -t4, t3, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (28->10), mult. (30->11), div. (0->0), fcn. (34->5), ass. (0->10)
	t11 = pkin(1) + r_i_i_C(3) + pkin(7) + qJ(3);
	t3 = pkin(9) + qJ(4);
	t1 = sin(t3);
	t2 = cos(t3);
	t10 = r_i_i_C(1) * t2 - r_i_i_C(2) * t1;
	t9 = -t1 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t8 = pkin(3) * sin(pkin(9)) + qJ(2) - t9;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t4 = [-t11 * t6 + t8 * t7, t6, t7, t10 * t6, 0, 0; t11 * t7 + t8 * t6, -t7, t6, -t10 * t7, 0, 0; 0, 0, 0, t9, 0, 0;];
	Ja_transl = t4;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (64->17), mult. (71->17), div. (0->0), fcn. (82->7), ass. (0->14)
	t3 = pkin(9) + qJ(4);
	t1 = sin(t3);
	t4 = sin(pkin(10));
	t6 = cos(pkin(10));
	t12 = r_i_i_C(1) * t6 - r_i_i_C(2) * t4 + pkin(4);
	t13 = r_i_i_C(3) + qJ(5);
	t2 = cos(t3);
	t15 = t13 * t1 + t12 * t2;
	t14 = -t12 * t1 + t13 * t2;
	t11 = t4 * r_i_i_C(1) + t6 * r_i_i_C(2) + pkin(1) + pkin(7) + qJ(3);
	t10 = sin(pkin(9)) * pkin(3) + qJ(2) - t14;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t5 = [t10 * t9 - t11 * t8, t8, t9, t15 * t8, -t8 * t2, 0; t10 * t8 + t11 * t9, -t9, t8, -t15 * t9, t9 * t2, 0; 0, 0, 0, t14, t1, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (113->30), mult. (100->37), div. (0->0), fcn. (115->9), ass. (0->24)
	t21 = r_i_i_C(3) + pkin(8) + qJ(5);
	t11 = pkin(9) + qJ(4);
	t9 = cos(t11);
	t27 = t21 * t9;
	t5 = cos(pkin(10)) * pkin(5) + pkin(4);
	t10 = pkin(10) + qJ(6);
	t6 = sin(t10);
	t8 = cos(t10);
	t19 = r_i_i_C(1) * t8 - r_i_i_C(2) * t6 + t5;
	t7 = sin(t11);
	t26 = t19 * t9 + t21 * t7;
	t16 = sin(qJ(1));
	t25 = t16 * t6;
	t24 = t16 * t8;
	t17 = cos(qJ(1));
	t23 = t17 * t6;
	t22 = t17 * t8;
	t20 = pkin(5) * sin(pkin(10)) + pkin(1) + pkin(7) + qJ(3);
	t18 = pkin(3) * sin(pkin(9)) - t27 + t7 * t5 + qJ(2);
	t4 = t7 * t22 - t25;
	t3 = t7 * t23 + t24;
	t2 = t7 * t24 + t23;
	t1 = -t7 * t25 + t22;
	t12 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t20 * t16 + t18 * t17, t16, t17, t26 * t16, -t16 * t9, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t18 * t16 + t20 * t17, -t17, t16, -t26 * t17, t17 * t9, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 0, 0, -t19 * t7 + t27, t7, (-r_i_i_C(1) * t6 - r_i_i_C(2) * t8) * t9;];
	Ja_transl = t12;
else
	Ja_transl=NaN(3,6);
end