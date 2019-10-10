% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:29
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPPRR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
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
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
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
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (27->14), mult. (44->14), div. (0->0), fcn. (60->6), ass. (0->10)
	t13 = pkin(1) + pkin(2);
	t12 = r_i_i_C(3) + qJ(4);
	t11 = cos(pkin(9));
	t10 = sin(pkin(9));
	t9 = r_i_i_C(1) * cos(pkin(10)) - r_i_i_C(2) * sin(pkin(10)) + pkin(3);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t2 = t10 * t8 - t11 * t7;
	t1 = -t10 * t7 - t11 * t8;
	t3 = [t8 * qJ(2) + t1 * t12 - t13 * t7 + t2 * t9, t7, 0, t2, 0, 0; t7 * qJ(2) - t1 * t9 + t12 * t2 + t13 * t8, -t8, 0, -t1, 0, 0; 0, 0, -1, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (48->17), mult. (64->19), div. (0->0), fcn. (84->7), ass. (0->15)
	t19 = pkin(1) + pkin(2);
	t18 = r_i_i_C(3) + pkin(7) + qJ(4);
	t17 = cos(qJ(1));
	t16 = sin(qJ(1));
	t15 = cos(pkin(9));
	t14 = sin(pkin(9));
	t9 = pkin(10) + qJ(5);
	t7 = sin(t9);
	t8 = cos(t9);
	t13 = -t8 * r_i_i_C(1) + t7 * r_i_i_C(2);
	t12 = r_i_i_C(1) * t7 + r_i_i_C(2) * t8;
	t11 = -t13 + cos(pkin(10)) * pkin(4) + pkin(3);
	t2 = t17 * t14 - t16 * t15;
	t1 = -t16 * t14 - t17 * t15;
	t3 = [t17 * qJ(2) + t18 * t1 + t11 * t2 - t19 * t16, t16, 0, t2, t12 * t1, 0; t16 * qJ(2) - t11 * t1 + t19 * t17 + t18 * t2, -t17, 0, -t1, t12 * t2, 0; 0, 0, -1, 0, t13, 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (115->32), mult. (167->41), div. (0->0), fcn. (217->9), ass. (0->25)
	t11 = pkin(10) + qJ(5);
	t10 = cos(t11);
	t13 = sin(qJ(6));
	t14 = cos(qJ(6));
	t16 = t14 * r_i_i_C(1) - t13 * r_i_i_C(2) + pkin(5);
	t25 = pkin(8) + r_i_i_C(3);
	t9 = sin(t11);
	t18 = t25 * t9;
	t28 = -t16 * t10 - t18;
	t27 = pkin(1) + pkin(2);
	t24 = cos(qJ(1));
	t23 = sin(qJ(1));
	t22 = t10 * t13;
	t21 = t10 * t14;
	t20 = cos(pkin(9));
	t19 = sin(pkin(9));
	t17 = t13 * r_i_i_C(1) + t14 * r_i_i_C(2);
	t15 = -t25 * t10 + t16 * t9;
	t12 = -pkin(7) - qJ(4);
	t8 = cos(pkin(10)) * pkin(4) + pkin(3);
	t4 = t24 * t19 - t23 * t20;
	t3 = -t23 * t19 - t24 * t20;
	t2 = t4 * t13 - t3 * t21;
	t1 = t14 * t4 + t3 * t22;
	t5 = [t24 * qJ(2) + (-t12 + t17) * t3 + (t8 - t28) * t4 - t27 * t23, t23, 0, t4, t15 * t3, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t23 * qJ(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t4 * t12 + (-t10 * pkin(5) - t18 - t8) * t3 + t27 * t24, -t24, 0, -t3, t15 * t4, (-t14 * t3 + t4 * t22) * r_i_i_C(1) + (t3 * t13 + t4 * t21) * r_i_i_C(2); 0, 0, -1, 0, t28, t17 * t9;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end