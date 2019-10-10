% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPPRR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
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
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
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
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (23->13), mult. (34->12), div. (0->0), fcn. (48->4), ass. (0->10)
	t12 = pkin(1) + pkin(2);
	t11 = pkin(3) - r_i_i_C(2);
	t10 = r_i_i_C(3) + qJ(4);
	t9 = cos(pkin(9));
	t8 = sin(pkin(9));
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t2 = -t6 * t9 + t7 * t8;
	t1 = -t6 * t8 - t7 * t9;
	t3 = [t7 * qJ(2) + t1 * t10 + t11 * t2 - t12 * t6, t6, 0, t2, 0, 0; t6 * qJ(2) - t1 * t11 + t10 * t2 + t12 * t7, -t7, 0, -t1, 0, 0; 0, 0, -1, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (38->16), mult. (68->18), div. (0->0), fcn. (90->6), ass. (0->14)
	t16 = pkin(1) + pkin(2);
	t15 = cos(qJ(1));
	t14 = sin(qJ(1));
	t13 = cos(pkin(9));
	t12 = sin(pkin(9));
	t11 = pkin(3) + pkin(7) + r_i_i_C(3);
	t6 = sin(qJ(5));
	t7 = cos(qJ(5));
	t10 = r_i_i_C(1) * t7 - r_i_i_C(2) * t6;
	t9 = t6 * r_i_i_C(1) + t7 * r_i_i_C(2);
	t8 = qJ(4) + t9;
	t2 = t15 * t12 - t14 * t13;
	t1 = -t14 * t12 - t15 * t13;
	t3 = [t15 * qJ(2) + t8 * t1 + t11 * t2 - t16 * t14, t14, 0, t2, t10 * t2, 0; t14 * qJ(2) - t11 * t1 + t16 * t15 + t8 * t2, -t15, 0, -t1, -t10 * t1, 0; 0, 0, -1, 0, t9, 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (84->31), mult. (171->40), div. (0->0), fcn. (223->8), ass. (0->23)
	t25 = pkin(1) + pkin(2);
	t11 = cos(qJ(5));
	t10 = cos(qJ(6));
	t8 = sin(qJ(6));
	t13 = t10 * r_i_i_C(1) - t8 * r_i_i_C(2) + pkin(5);
	t21 = pkin(8) + r_i_i_C(3);
	t9 = sin(qJ(5));
	t24 = t13 * t11 + t21 * t9;
	t23 = t8 * t9;
	t22 = -pkin(3) - pkin(7);
	t20 = t10 * t9;
	t19 = cos(qJ(1));
	t18 = sin(qJ(1));
	t17 = cos(pkin(9));
	t16 = sin(pkin(9));
	t15 = t21 * t11;
	t14 = t8 * r_i_i_C(1) + t10 * r_i_i_C(2);
	t12 = t13 * t9 - t15;
	t4 = t19 * t16 - t18 * t17;
	t3 = -t18 * t16 - t19 * t17;
	t2 = t4 * t20 - t3 * t8;
	t1 = -t10 * t3 - t4 * t23;
	t5 = [t19 * qJ(2) + (t14 - t22) * t4 + (qJ(4) + t12) * t3 - t25 * t18, t18, 0, t4, t24 * t4, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t18 * qJ(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t22 * t3 + (t9 * pkin(5) + qJ(4) - t15) * t4 + t25 * t19, -t19, 0, -t3, -t24 * t3, (-t10 * t4 + t3 * t23) * r_i_i_C(1) + (t3 * t20 + t4 * t8) * r_i_i_C(2); 0, 0, -1, 0, t12, t14 * t11;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end