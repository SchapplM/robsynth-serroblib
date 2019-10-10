% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:32
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPPRR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
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
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
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
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (9->5), mult. (10->4), div. (0->0), fcn. (14->2), ass. (0->5)
	t4 = r_i_i_C(1) + qJ(2);
	t3 = pkin(1) + r_i_i_C(3) + qJ(3);
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t5 = [-t3 * t1 + t4 * t2, t1, t2, 0, 0, 0; t4 * t1 + t3 * t2, -t2, t1, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (15->11), mult. (20->12), div. (0->0), fcn. (28->4), ass. (0->9)
	t8 = pkin(1) + qJ(3);
	t7 = pkin(3) + qJ(2);
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t4 = cos(pkin(9));
	t3 = sin(pkin(9));
	t2 = t6 * t3 + t5 * t4;
	t1 = -t5 * t3 + t6 * t4;
	t9 = [t1 * r_i_i_C(1) - t2 * r_i_i_C(2) - t8 * t5 + t7 * t6, t5, t6, 0, 0, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * t5 + t8 * t6, -t6, t5, 0, 0, 0; 0, 0, 0, 1, 0, 0;];
	Ja_transl = t9;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (34->16), mult. (60->18), div. (0->0), fcn. (78->6), ass. (0->15)
	t15 = pkin(7) + r_i_i_C(3);
	t14 = pkin(1) + qJ(3);
	t13 = pkin(3) + qJ(2);
	t12 = sin(pkin(9));
	t5 = sin(qJ(5));
	t7 = cos(qJ(5));
	t11 = t7 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t10 = r_i_i_C(1) * t5 + r_i_i_C(2) * t7;
	t9 = pkin(4) + t11;
	t8 = cos(qJ(1));
	t6 = sin(qJ(1));
	t4 = cos(pkin(9));
	t2 = t8 * t12 + t6 * t4;
	t1 = -t6 * t12 + t8 * t4;
	t3 = [t9 * t1 + t13 * t8 - t14 * t6 + t15 * t2, t6, t8, 0, -t10 * t2, 0; -t15 * t1 + t13 * t6 + t14 * t8 + t9 * t2, -t8, t6, 0, t10 * t1, 0; 0, 0, 0, 1, t11, 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (80->30), mult. (163->40), div. (0->0), fcn. (211->8), ass. (0->23)
	t23 = pkin(1) + qJ(3);
	t22 = pkin(3) + qJ(2);
	t10 = cos(qJ(5));
	t7 = sin(qJ(6));
	t9 = cos(qJ(6));
	t13 = t9 * r_i_i_C(1) - t7 * r_i_i_C(2) + pkin(5);
	t20 = pkin(8) + r_i_i_C(3);
	t8 = sin(qJ(5));
	t21 = t20 * t10 - t13 * t8;
	t19 = t10 * t7;
	t18 = t10 * t9;
	t17 = sin(qJ(1));
	t16 = sin(pkin(9));
	t15 = t20 * t8;
	t14 = -t7 * r_i_i_C(1) - t9 * r_i_i_C(2);
	t12 = t13 * t10 + t15;
	t11 = cos(qJ(1));
	t6 = cos(pkin(9));
	t4 = t11 * t16 + t17 * t6;
	t3 = t11 * t6 - t17 * t16;
	t2 = t4 * t18 - t3 * t7;
	t1 = -t4 * t19 - t3 * t9;
	t5 = [t22 * t11 + (pkin(7) - t14) * t4 + (pkin(4) + t12) * t3 - t23 * t17, t17, t11, 0, t21 * t4, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; -t3 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t23 * t11 + (pkin(5) * t10 + pkin(4) + t15) * t4 + t22 * t17, -t11, t17, 0, -t21 * t3, (t3 * t19 - t4 * t9) * r_i_i_C(1) + (t3 * t18 + t4 * t7) * r_i_i_C(2); 0, 0, 0, 1, t12, t14 * t8;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end