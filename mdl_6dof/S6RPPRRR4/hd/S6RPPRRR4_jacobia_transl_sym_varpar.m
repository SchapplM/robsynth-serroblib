% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:06
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
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
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
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
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (13->10), mult. (18->12), div. (0->0), fcn. (24->4), ass. (0->8)
	t7 = pkin(1) + pkin(2);
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t4 = cos(pkin(10));
	t3 = sin(pkin(10));
	t2 = t3 * t5 + t4 * t6;
	t1 = t3 * t6 - t4 * t5;
	t8 = [t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + t6 * qJ(2) - t5 * t7, t5, 0, 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + t5 * qJ(2) + t6 * t7, -t6, 0, 0, 0, 0; 0, 0, -1, 0, 0, 0;];
	Ja_transl = t8;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (32->14), mult. (58->18), div. (0->0), fcn. (74->6), ass. (0->14)
	t16 = pkin(1) + pkin(2);
	t15 = pkin(7) + r_i_i_C(3);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t12 = cos(pkin(10));
	t11 = sin(pkin(10));
	t6 = sin(qJ(4));
	t7 = cos(qJ(4));
	t10 = -t7 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t9 = r_i_i_C(1) * t6 + r_i_i_C(2) * t7;
	t8 = pkin(3) - t10;
	t2 = t11 * t14 - t12 * t13;
	t1 = -t13 * t11 - t14 * t12;
	t3 = [t14 * qJ(2) + t15 * t1 - t16 * t13 + t8 * t2, t13, 0, t9 * t1, 0, 0; t13 * qJ(2) - t8 * t1 + t16 * t14 + t15 * t2, -t14, 0, t9 * t2, 0, 0; 0, 0, -1, t10, 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (78->28), mult. (161->40), div. (0->0), fcn. (207->8), ass. (0->22)
	t11 = cos(qJ(4));
	t10 = cos(qJ(5));
	t8 = sin(qJ(5));
	t13 = r_i_i_C(1) * t10 - r_i_i_C(2) * t8 + pkin(4);
	t22 = pkin(8) + r_i_i_C(3);
	t9 = sin(qJ(4));
	t15 = t22 * t9;
	t25 = -t13 * t11 - t15;
	t24 = pkin(1) + pkin(2);
	t21 = t11 * t8;
	t20 = cos(qJ(1));
	t19 = sin(qJ(1));
	t18 = t10 * t11;
	t17 = cos(pkin(10));
	t16 = sin(pkin(10));
	t14 = t8 * r_i_i_C(1) + t10 * r_i_i_C(2);
	t12 = -t22 * t11 + t13 * t9;
	t4 = t20 * t16 - t19 * t17;
	t3 = -t19 * t16 - t20 * t17;
	t2 = -t3 * t18 + t4 * t8;
	t1 = t4 * t10 + t3 * t21;
	t5 = [t20 * qJ(2) + (pkin(7) + t14) * t3 + (pkin(3) - t25) * t4 - t24 * t19, t19, 0, t12 * t3, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t19 * qJ(2) + t4 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + (-pkin(4) * t11 - pkin(3) - t15) * t3 + t24 * t20, -t20, 0, t12 * t4, (-t3 * t10 + t4 * t21) * r_i_i_C(1) + (t4 * t18 + t3 * t8) * r_i_i_C(2), 0; 0, 0, -1, t25, t14 * t9, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (150->36), mult. (230->50), div. (0->0), fcn. (294->10), ass. (0->32)
	t21 = cos(qJ(4));
	t20 = cos(qJ(5));
	t14 = t20 * pkin(5) + pkin(4);
	t17 = qJ(5) + qJ(6);
	t15 = sin(t17);
	t16 = cos(t17);
	t24 = r_i_i_C(1) * t16 - r_i_i_C(2) * t15 + t14;
	t19 = sin(qJ(4));
	t35 = r_i_i_C(3) + pkin(9) + pkin(8);
	t25 = t35 * t19;
	t44 = -t24 * t21 - t25;
	t43 = t15 * r_i_i_C(1) + t16 * r_i_i_C(2);
	t42 = pkin(1) + pkin(2);
	t30 = t16 * t21;
	t31 = t15 * t21;
	t27 = sin(pkin(10));
	t28 = cos(pkin(10));
	t32 = sin(qJ(1));
	t33 = cos(qJ(1));
	t7 = -t32 * t27 - t33 * t28;
	t8 = t33 * t27 - t32 * t28;
	t40 = (-t7 * t16 + t8 * t31) * r_i_i_C(1) + (t7 * t15 + t8 * t30) * r_i_i_C(2);
	t5 = t8 * t16 + t7 * t31;
	t6 = t8 * t15 - t7 * t30;
	t39 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t18 = sin(qJ(5));
	t38 = pkin(5) * t18;
	t34 = t43 * t19;
	t29 = t18 * t21;
	t26 = pkin(7) + t38;
	t23 = t24 * t19 - t35 * t21;
	t1 = [t33 * qJ(2) + (t26 + t43) * t7 + (pkin(3) - t44) * t8 - t42 * t32, t32, 0, t23 * t7, (t20 * t8 + t7 * t29) * pkin(5) + t39, t39; t32 * qJ(2) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t26 * t8 + (-t21 * t14 - pkin(3) - t25) * t7 + t42 * t33, -t33, 0, t23 * t8, (-t20 * t7 + t8 * t29) * pkin(5) + t40, t40; 0, 0, -1, t44, t19 * t38 + t34, t34;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end