% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:53
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
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
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (9->6), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(10)) - r_i_i_C(2) * sin(pkin(10)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t5 * t3 + t6 * t4, t3, 0, 0, 0, 0; t6 * t3 + t5 * t4, -t4, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (26->9), mult. (26->11), div. (0->0), fcn. (28->5), ass. (0->10)
	t11 = r_i_i_C(3) + pkin(7) + qJ(2);
	t4 = pkin(10) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t10 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t9 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t8 = cos(pkin(10)) * pkin(2) + pkin(1) + t10;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [t11 * t7 - t8 * t6, t6, t9 * t7, 0, 0, 0; t11 * t6 + t8 * t7, -t7, t9 * t6, 0, 0, 0; 0, 0, t10, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (50->12), mult. (46->13), div. (0->0), fcn. (51->5), ass. (0->12)
	t10 = r_i_i_C(3) + qJ(4);
	t12 = pkin(3) + r_i_i_C(1);
	t4 = pkin(10) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t9 = t10 * t2 + t12 * t3;
	t13 = cos(pkin(10)) * pkin(2) + pkin(1) + t9;
	t11 = r_i_i_C(2) + pkin(7) + qJ(2);
	t8 = t10 * t3 - t12 * t2;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [t11 * t7 - t13 * t6, t6, t8 * t7, t7 * t2, 0, 0; t11 * t6 + t13 * t7, -t7, t8 * t6, t6 * t2, 0, 0; 0, 0, t9, -t3, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (103->24), mult. (111->32), div. (0->0), fcn. (132->7), ass. (0->24)
	t11 = pkin(10) + qJ(3);
	t10 = cos(t11);
	t28 = pkin(3) + pkin(4);
	t9 = sin(t11);
	t17 = t9 * qJ(4) + t28 * t10;
	t29 = cos(pkin(10)) * pkin(2) + pkin(1) + t17;
	t16 = cos(qJ(1));
	t27 = t16 * t9;
	t13 = sin(qJ(5));
	t26 = t10 * t13;
	t24 = -pkin(8) - r_i_i_C(3) + pkin(7) + qJ(2);
	t14 = sin(qJ(1));
	t15 = cos(qJ(5));
	t19 = -t15 * t9 + t26;
	t1 = t19 * t14;
	t20 = t10 * t15 + t13 * t9;
	t2 = t20 * t14;
	t23 = -t1 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t3 = -t15 * t27 + t16 * t26;
	t4 = t20 * t16;
	t22 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t21 = -r_i_i_C(1) * t20 + r_i_i_C(2) * t19;
	t18 = qJ(4) * t10 - t28 * t9;
	t5 = [-t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t29 * t14 + t24 * t16, t14, t18 * t16 - t22, t27, t22, 0; t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t24 * t14 + t29 * t16, -t16, t18 * t14 - t23, t14 * t9, t23, 0; 0, 0, t17 - t21, -t10, t21, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (226->38), mult. (273->54), div. (0->0), fcn. (332->9), ass. (0->31)
	t17 = pkin(10) + qJ(3);
	t15 = sin(t17);
	t16 = cos(t17);
	t36 = pkin(3) + pkin(4);
	t23 = t15 * qJ(4) + t36 * t16;
	t40 = cos(pkin(10)) * pkin(2) + pkin(1) + t23;
	t32 = sin(qJ(5));
	t33 = cos(qJ(5));
	t8 = t15 * t33 - t16 * t32;
	t19 = sin(qJ(6));
	t21 = cos(qJ(6));
	t25 = t21 * r_i_i_C(1) - t19 * r_i_i_C(2) + pkin(5);
	t20 = sin(qJ(1));
	t3 = t8 * t20;
	t35 = -r_i_i_C(3) - pkin(9);
	t7 = t15 * t32 + t16 * t33;
	t4 = t7 * t20;
	t39 = t25 * t3 - t35 * t4;
	t22 = cos(qJ(1));
	t27 = t22 * t32;
	t28 = t22 * t33;
	t5 = -t15 * t27 - t16 * t28;
	t6 = -t15 * t28 + t16 * t27;
	t38 = -t25 * t6 + t35 * t5;
	t37 = -t25 * t7 - t35 * t8;
	t34 = -pkin(8) + pkin(7) + qJ(2);
	t26 = -t19 * r_i_i_C(1) - t21 * r_i_i_C(2);
	t24 = qJ(4) * t16 - t36 * t15;
	t2 = -t20 * t19 - t5 * t21;
	t1 = t5 * t19 - t20 * t21;
	t9 = [-t35 * t3 - t25 * t4 + (t26 + t34) * t22 - t40 * t20, t20, t24 * t22 - t38, t22 * t15, t38, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); -t5 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t34 * t20 + t40 * t22 - t35 * t6, -t22, t24 * t20 - t39, t20 * t15, t39, (-t4 * t19 + t22 * t21) * r_i_i_C(1) + (-t22 * t19 - t4 * t21) * r_i_i_C(2); 0, 0, t23 - t37, -t16, t37, t26 * t8;];
	Ja_transl = t9;
else
	Ja_transl=NaN(3,6);
end