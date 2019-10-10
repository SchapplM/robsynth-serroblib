% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRPR6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:42
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRPR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
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
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
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
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (9->5), mult. (10->4), div. (0->0), fcn. (14->2), ass. (0->5)
	t4 = r_i_i_C(2) + qJ(2);
	t3 = pkin(1) + r_i_i_C(3) + qJ(3);
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t5 = [-t3 * t1 + t4 * t2, t1, t2, 0, 0, 0; t4 * t1 + t3 * t2, -t2, t1, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (16->7), mult. (28->10), div. (0->0), fcn. (32->4), ass. (0->9)
	t8 = -pkin(7) - r_i_i_C(3) + qJ(2);
	t1 = sin(qJ(4));
	t3 = cos(qJ(4));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -t1 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t5 = pkin(1) + qJ(3) - t6;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t5 * t2 + t8 * t4, t2, t4, t7 * t4, 0, 0; t8 * t2 + t5 * t4, -t4, t2, t7 * t2, 0, 0; 0, 0, 0, t6, 0, 0;];
	Ja_transl = t9;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (28->11), mult. (48->12), div. (0->0), fcn. (55->4), ass. (0->11)
	t1 = sin(qJ(4));
	t10 = pkin(4) - r_i_i_C(2);
	t3 = cos(qJ(4));
	t8 = r_i_i_C(3) + qJ(5);
	t5 = -t10 * t1 + t8 * t3;
	t11 = -pkin(1) - qJ(3) + t5;
	t7 = -pkin(7) - r_i_i_C(1) + qJ(2);
	t6 = t8 * t1 + t10 * t3;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [t11 * t2 + t7 * t4, t2, t4, t6 * t4, -t4 * t3, 0; -t11 * t4 + t7 * t2, -t4, t2, t6 * t2, -t2 * t3, 0; 0, 0, 0, t5, t1, 0;];
	Ja_transl = t9;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (49->25), mult. (100->34), div. (0->0), fcn. (115->6), ass. (0->19)
	t16 = pkin(4) + pkin(8) + r_i_i_C(3);
	t6 = sin(qJ(4));
	t14 = t16 * t6;
	t9 = cos(qJ(4));
	t19 = -t9 * qJ(5) + pkin(1) + qJ(3) + t14;
	t7 = sin(qJ(1));
	t18 = t7 * t9;
	t10 = cos(qJ(1));
	t17 = t10 * t9;
	t15 = -pkin(5) - pkin(7) + qJ(2);
	t5 = sin(qJ(6));
	t8 = cos(qJ(6));
	t12 = r_i_i_C(1) * t5 + r_i_i_C(2) * t8 + qJ(5);
	t11 = t12 * t6 + t16 * t9;
	t4 = -t10 * t8 + t5 * t18;
	t3 = t10 * t5 + t8 * t18;
	t2 = t5 * t17 + t7 * t8;
	t1 = -t8 * t17 + t7 * t5;
	t13 = [t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t15 * t10 - t19 * t7, t7, t10, t11 * t10, -t17, r_i_i_C(1) * t1 + r_i_i_C(2) * t2; -t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t19 * t10 + t15 * t7, -t10, t7, t11 * t7, -t18, -r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 0, 0, t12 * t9 - t14, t6, (r_i_i_C(1) * t8 - r_i_i_C(2) * t5) * t6;];
	Ja_transl = t13;
else
	Ja_transl=NaN(3,6);
end