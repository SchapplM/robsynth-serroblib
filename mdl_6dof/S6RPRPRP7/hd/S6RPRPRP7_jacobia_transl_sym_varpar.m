% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRP7_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP7_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP7_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:23
	% EndTime: 2019-10-10 00:39:23
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:23
	% EndTime: 2019-10-10 00:39:23
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
	% StartTime: 2019-10-10 00:39:23
	% EndTime: 2019-10-10 00:39:23
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-10 00:39:23
	% EndTime: 2019-10-10 00:39:23
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (14->8), mult. (26->10), div. (0->0), fcn. (28->4), ass. (0->9)
	t8 = pkin(1) + pkin(7) + r_i_i_C(3);
	t1 = sin(qJ(3));
	t3 = cos(qJ(3));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -t1 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t5 = qJ(2) - t6;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t8 * t2 + t5 * t4, t2, t7 * t2, 0, 0, 0; t5 * t2 + t8 * t4, -t4, -t7 * t4, 0, 0, 0; 0, 0, t6, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:23
	% EndTime: 2019-10-10 00:39:23
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (31->11), mult. (35->12), div. (0->0), fcn. (39->6), ass. (0->10)
	t12 = pkin(1) + r_i_i_C(3) + qJ(4) + pkin(7);
	t3 = qJ(3) + pkin(9);
	t1 = sin(t3);
	t2 = cos(t3);
	t11 = -t1 * r_i_i_C(1) - t2 * r_i_i_C(2) - sin(qJ(3)) * pkin(3);
	t10 = cos(qJ(3)) * pkin(3) + r_i_i_C(1) * t2 - r_i_i_C(2) * t1;
	t9 = qJ(2) - t11;
	t8 = cos(qJ(1));
	t6 = sin(qJ(1));
	t4 = [-t12 * t6 + t8 * t9, t6, t10 * t6, t8, 0, 0; t12 * t8 + t6 * t9, -t8, -t10 * t8, t6, 0, 0; 0, 0, t11, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:23
	% EndTime: 2019-10-10 00:39:23
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (76->26), mult. (94->34), div. (0->0), fcn. (106->8), ass. (0->22)
	t23 = pkin(8) + r_i_i_C(3);
	t7 = qJ(3) + pkin(9);
	t6 = cos(t7);
	t26 = t23 * t6 - sin(qJ(3)) * pkin(3);
	t12 = cos(qJ(5));
	t9 = sin(qJ(5));
	t16 = r_i_i_C(1) * t12 - r_i_i_C(2) * t9 + pkin(4);
	t5 = sin(t7);
	t25 = t16 * t6 + t23 * t5 + cos(qJ(3)) * pkin(3);
	t24 = pkin(1) + qJ(4) + pkin(7);
	t11 = sin(qJ(1));
	t20 = t11 * t9;
	t14 = cos(qJ(1));
	t19 = t14 * t9;
	t18 = t11 * t12;
	t17 = t12 * t14;
	t15 = t5 * pkin(4) + qJ(2) - t26;
	t4 = t5 * t17 - t20;
	t3 = t5 * t19 + t18;
	t2 = t5 * t18 + t19;
	t1 = -t5 * t20 + t17;
	t8 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t24 * t11 + t15 * t14, t11, t25 * t11, t14, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * t11 + t24 * t14, -t14, -t25 * t14, t11, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0; 0, 0, -t16 * t5 + t26, 0, (-r_i_i_C(1) * t9 - r_i_i_C(2) * t12) * t6, 0;];
	Ja_transl = t8;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:23
	% EndTime: 2019-10-10 00:39:23
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (100->31), mult. (117->38), div. (0->0), fcn. (132->8), ass. (0->24)
	t29 = pkin(5) + r_i_i_C(1);
	t26 = r_i_i_C(3) + qJ(6) + pkin(8);
	t8 = qJ(3) + pkin(9);
	t7 = cos(t8);
	t28 = t26 * t7 - sin(qJ(3)) * pkin(3);
	t11 = sin(qJ(5));
	t14 = cos(qJ(5));
	t5 = pkin(5) * t14 + pkin(4);
	t18 = r_i_i_C(1) * t14 - r_i_i_C(2) * t11 + t5;
	t6 = sin(t8);
	t27 = t18 * t7 + t26 * t6 + cos(qJ(3)) * pkin(3);
	t16 = cos(qJ(1));
	t23 = t11 * t16;
	t13 = sin(qJ(1));
	t22 = t13 * t11;
	t21 = t13 * t14;
	t20 = t14 * t16;
	t19 = pkin(5) * t11 + pkin(1) + pkin(7) + qJ(4);
	t3 = t6 * t23 + t21;
	t1 = -t6 * t22 + t20;
	t17 = t5 * t6 + qJ(2) - t28;
	t4 = t6 * t20 - t22;
	t2 = t6 * t21 + t23;
	t9 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t19 * t13 + t17 * t16, t13, t27 * t13, t16, -t2 * r_i_i_C(2) + t29 * t1, -t13 * t7; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t17 * t13 + t19 * t16, -t16, -t27 * t16, t13, t4 * r_i_i_C(2) + t29 * t3, t16 * t7; 0, 0, -t18 * t6 + t28, 0, (-r_i_i_C(2) * t14 - t29 * t11) * t7, t6;];
	Ja_transl = t9;
else
	Ja_transl=NaN(3,6);
end