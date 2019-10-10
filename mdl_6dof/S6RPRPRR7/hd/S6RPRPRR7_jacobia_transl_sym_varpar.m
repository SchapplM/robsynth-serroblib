% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:56
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR7_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR7_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR7_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
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
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
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
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
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
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (31->11), mult. (35->12), div. (0->0), fcn. (39->6), ass. (0->10)
	t12 = pkin(1) + r_i_i_C(3) + qJ(4) + pkin(7);
	t3 = qJ(3) + pkin(10);
	t1 = sin(t3);
	t2 = cos(t3);
	t11 = -t1 * r_i_i_C(1) - t2 * r_i_i_C(2) - sin(qJ(3)) * pkin(3);
	t10 = cos(qJ(3)) * pkin(3) + r_i_i_C(1) * t2 - r_i_i_C(2) * t1;
	t9 = qJ(2) - t11;
	t8 = cos(qJ(1));
	t6 = sin(qJ(1));
	t4 = [-t12 * t6 + t9 * t8, t6, t10 * t6, t8, 0, 0; t12 * t8 + t9 * t6, -t8, -t10 * t8, t6, 0, 0; 0, 0, t11, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (68->18), mult. (50->18), div. (0->0), fcn. (54->8), ass. (0->16)
	t10 = qJ(3) + pkin(10);
	t8 = qJ(5) + t10;
	t5 = sin(t8);
	t6 = cos(t8);
	t15 = t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t19 = -sin(qJ(3)) * pkin(3) - pkin(4) * sin(t10) - t15;
	t18 = r_i_i_C(1) * t6;
	t17 = r_i_i_C(2) * t5;
	t16 = pkin(1) + r_i_i_C(3) + pkin(8) + qJ(4) + pkin(7);
	t14 = qJ(2) - t19;
	t13 = cos(qJ(1));
	t12 = sin(qJ(1));
	t4 = t13 * t17;
	t3 = t12 * t18;
	t2 = pkin(4) * cos(t10) + cos(qJ(3)) * pkin(3);
	t1 = [-t12 * t16 + t14 * t13, t12, t3 + (t2 - t17) * t12, t13, -t12 * t17 + t3, 0; t14 * t12 + t13 * t16, -t13, t4 + (-t2 - t18) * t13, t12, -t13 * t18 + t4, 0; 0, 0, t19, 0, -t15, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (161->36), mult. (130->44), div. (0->0), fcn. (142->10), ass. (0->31)
	t20 = qJ(3) + pkin(10);
	t18 = qJ(5) + t20;
	t15 = sin(t18);
	t40 = pkin(9) + r_i_i_C(3);
	t42 = t40 * t15;
	t16 = cos(t18);
	t41 = t40 * t16;
	t21 = sin(qJ(6));
	t38 = r_i_i_C(2) * t21;
	t37 = pkin(1) + pkin(8) + qJ(4) + pkin(7);
	t25 = cos(qJ(1));
	t35 = t21 * t25;
	t23 = sin(qJ(1));
	t34 = t23 * t21;
	t24 = cos(qJ(6));
	t33 = t23 * t24;
	t32 = t24 * t25;
	t31 = t16 * t38;
	t30 = -r_i_i_C(1) * t24 - pkin(5);
	t29 = t23 * t42 + (pkin(5) * t23 + r_i_i_C(1) * t33) * t16;
	t7 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t20);
	t28 = t41 + (t30 + t38) * t15;
	t27 = pkin(5) * t15 + qJ(2) - t41 + t7;
	t26 = t30 * t16 - t42;
	t8 = pkin(4) * cos(t20) + cos(qJ(3)) * pkin(3);
	t6 = t25 * t31;
	t4 = t15 * t32 - t34;
	t3 = t15 * t35 + t33;
	t2 = t15 * t33 + t35;
	t1 = -t15 * t34 + t32;
	t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t37 * t23 + t27 * t25, t23, (t8 - t31) * t23 + t29, t25, -t23 * t31 + t29, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t27 * t23 + t37 * t25, -t25, t6 + (-t8 + t26) * t25, t23, t26 * t25 + t6, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 0, t28 - t7, 0, t28, (-r_i_i_C(1) * t21 - r_i_i_C(2) * t24) * t16;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end