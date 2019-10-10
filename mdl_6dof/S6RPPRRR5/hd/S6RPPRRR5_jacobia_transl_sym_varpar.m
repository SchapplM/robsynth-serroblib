% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:08
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:18
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:18
	% EndTime: 2019-10-10 00:08:19
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
	% StartTime: 2019-10-10 00:08:18
	% EndTime: 2019-10-10 00:08:19
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
	% StartTime: 2019-10-10 00:08:18
	% EndTime: 2019-10-10 00:08:19
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
	% StartTime: 2019-10-10 00:08:18
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-10-10 00:08:18
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (42->13), mult. (47->16), div. (0->0), fcn. (51->6), ass. (0->15)
	t5 = qJ(4) + qJ(5);
	t4 = cos(t5);
	t17 = r_i_i_C(1) * t4;
	t3 = sin(t5);
	t16 = r_i_i_C(2) * t3;
	t15 = -r_i_i_C(3) + qJ(2) - pkin(8) - pkin(7);
	t14 = cos(qJ(4)) * pkin(4) - t16;
	t13 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t12 = -sin(qJ(4)) * pkin(4) + t13;
	t11 = pkin(1) + qJ(3) - t12;
	t9 = cos(qJ(1));
	t7 = sin(qJ(1));
	t2 = t9 * t17;
	t1 = t7 * t17;
	t6 = [-t11 * t7 + t15 * t9, t7, t9, t14 * t9 + t2, -t16 * t9 + t2, 0; t11 * t9 + t15 * t7, -t9, t7, t14 * t7 + t1, -t16 * t7 + t1, 0; 0, 0, 0, t12, t13, 0;];
	Ja_transl = t6;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:18
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (105->31), mult. (127->39), div. (0->0), fcn. (139->8), ass. (0->29)
	t22 = cos(qJ(6));
	t51 = r_i_i_C(1) * t22 + pkin(5);
	t18 = qJ(4) + qJ(5);
	t16 = sin(t18);
	t17 = cos(t18);
	t49 = pkin(9) + r_i_i_C(3);
	t50 = t16 * t49 + t51 * t17;
	t46 = t49 * t17;
	t43 = sin(qJ(4)) * pkin(4);
	t45 = -t16 * pkin(5) - pkin(1) - qJ(3) - t43 + t46;
	t19 = sin(qJ(6));
	t40 = r_i_i_C(2) * t19;
	t21 = sin(qJ(1));
	t37 = t19 * t21;
	t24 = cos(qJ(1));
	t36 = t19 * t24;
	t35 = t21 * t22;
	t34 = t22 * t24;
	t33 = qJ(2) - pkin(8) - pkin(7);
	t31 = t17 * t40;
	t30 = t50 * t21;
	t29 = t50 * t24;
	t28 = cos(qJ(4)) * pkin(4) - t31;
	t26 = t46 + (t40 - t51) * t16;
	t4 = t16 * t34 - t37;
	t3 = -t16 * t36 - t35;
	t2 = -t16 * t35 - t36;
	t1 = t16 * t37 - t34;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t45 * t21 + t33 * t24, t21, t24, t28 * t24 + t29, -t24 * t31 + t29, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t33 * t21 - t45 * t24, -t24, t21, t28 * t21 + t30, -t21 * t31 + t30, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, 0, 0, t26 - t43, t26, (-r_i_i_C(1) * t19 - r_i_i_C(2) * t22) * t17;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end