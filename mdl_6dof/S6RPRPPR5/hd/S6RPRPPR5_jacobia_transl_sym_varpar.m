% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:22
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPPR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (9->6), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(9)) - r_i_i_C(2) * sin(pkin(9)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t5 * t3 + t6 * t4, t3, 0, 0, 0, 0; t6 * t3 + t5 * t4, -t4, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (26->9), mult. (26->11), div. (0->0), fcn. (28->5), ass. (0->10)
	t11 = r_i_i_C(3) + pkin(7) + qJ(2);
	t4 = pkin(9) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t10 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t9 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t8 = cos(pkin(9)) * pkin(2) + pkin(1) + t10;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [t11 * t7 - t8 * t6, t6, t9 * t7, 0, 0, 0; t11 * t6 + t8 * t7, -t7, t9 * t6, 0, 0, 0; 0, 0, t10, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (50->12), mult. (46->13), div. (0->0), fcn. (51->5), ass. (0->12)
	t10 = r_i_i_C(3) + qJ(4);
	t12 = pkin(3) - r_i_i_C(2);
	t4 = pkin(9) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t9 = t10 * t2 + t12 * t3;
	t13 = cos(pkin(9)) * pkin(2) + pkin(1) + t9;
	t11 = r_i_i_C(1) + pkin(7) + qJ(2);
	t8 = t10 * t3 - t12 * t2;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [t11 * t7 - t13 * t6, t6, t8 * t7, t7 * t2, 0, 0; t11 * t6 + t13 * t7, -t7, t8 * t6, t6 * t2, 0, 0; 0, 0, t9, -t3, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (77->16), mult. (80->19), div. (0->0), fcn. (92->7), ass. (0->14)
	t5 = sin(pkin(10));
	t6 = cos(pkin(10));
	t13 = r_i_i_C(1) * t5 + r_i_i_C(2) * t6 + qJ(4);
	t14 = pkin(3) + r_i_i_C(3) + qJ(5);
	t4 = pkin(9) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t11 = t13 * t2 + t14 * t3;
	t15 = cos(pkin(9)) * pkin(2) + pkin(1) + t11;
	t12 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2) + pkin(4) + pkin(7) + qJ(2);
	t10 = t13 * t3 - t14 * t2;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [t12 * t9 - t15 * t8, t8, t10 * t9, t9 * t2, t9 * t3, 0; t12 * t8 + t15 * t9, -t9, t10 * t8, t8 * t2, t8 * t3, 0; 0, 0, t11, -t3, t2, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:18
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (131->29), mult. (116->39), div. (0->0), fcn. (132->9), ass. (0->24)
	t12 = pkin(9) + qJ(3);
	t10 = cos(t12);
	t22 = pkin(3) + r_i_i_C(3) + pkin(8) + qJ(5);
	t20 = t22 * t10;
	t21 = pkin(5) * sin(pkin(10)) + qJ(4);
	t8 = sin(t12);
	t29 = -t21 * t8 - cos(pkin(9)) * pkin(2) - pkin(1) - t20;
	t16 = sin(qJ(1));
	t27 = t16 * t8;
	t11 = pkin(10) + qJ(6);
	t9 = cos(t11);
	t26 = t16 * t9;
	t17 = cos(qJ(1));
	t25 = t17 * t8;
	t24 = t17 * t9;
	t23 = pkin(7) + qJ(2) + cos(pkin(10)) * pkin(5) + pkin(4);
	t7 = sin(t11);
	t19 = r_i_i_C(1) * t7 + r_i_i_C(2) * t9 + t21;
	t18 = t19 * t10 - t22 * t8;
	t4 = -t7 * t27 + t24;
	t3 = t17 * t7 + t8 * t26;
	t2 = t7 * t25 + t26;
	t1 = -t16 * t7 + t8 * t24;
	t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t29 * t16 + t23 * t17, t16, t18 * t17, t25, t17 * t10, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t23 * t16 - t29 * t17, -t17, t18 * t16, t27, t16 * t10, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 0, t19 * t8 + t20, -t10, t8, (-r_i_i_C(1) * t9 + r_i_i_C(2) * t7) * t10;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end