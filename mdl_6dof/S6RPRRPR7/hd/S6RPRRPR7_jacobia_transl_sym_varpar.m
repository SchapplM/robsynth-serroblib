% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:33
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR7_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR7_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR7_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:55
	% EndTime: 2019-10-10 01:33:55
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:55
	% EndTime: 2019-10-10 01:33:56
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
	% StartTime: 2019-10-10 01:33:55
	% EndTime: 2019-10-10 01:33:56
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
	% StartTime: 2019-10-10 01:33:55
	% EndTime: 2019-10-10 01:33:56
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
	% StartTime: 2019-10-10 01:33:55
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (40->14), mult. (45->16), div. (0->0), fcn. (47->6), ass. (0->15)
	t17 = cos(qJ(3)) * pkin(3);
	t5 = qJ(3) + qJ(4);
	t4 = cos(t5);
	t16 = r_i_i_C(1) * t4;
	t3 = sin(t5);
	t15 = r_i_i_C(2) * t3;
	t14 = pkin(1) + r_i_i_C(3) + pkin(8) + pkin(7);
	t13 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t12 = -sin(qJ(3)) * pkin(3) + t13;
	t11 = qJ(2) - t12;
	t9 = cos(qJ(1));
	t7 = sin(qJ(1));
	t2 = t9 * t15;
	t1 = t7 * t16;
	t6 = [t11 * t9 - t14 * t7, t7, t1 + (-t15 + t17) * t7, -t7 * t15 + t1, 0, 0; t11 * t7 + t14 * t9, -t9, t2 + (-t16 - t17) * t9, -t9 * t16 + t2, 0, 0; 0, 0, t12, t13, 0, 0;];
	Ja_transl = t6;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:55
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (74->19), mult. (55->18), div. (0->0), fcn. (59->8), ass. (0->17)
	t11 = qJ(3) + qJ(4);
	t7 = pkin(10) + t11;
	t5 = sin(t7);
	t6 = cos(t7);
	t15 = -t5 * r_i_i_C(1) - t6 * r_i_i_C(2) - pkin(4) * sin(t11);
	t24 = t15 - sin(qJ(3)) * pkin(3);
	t22 = pkin(4) * cos(t11);
	t21 = r_i_i_C(1) * t6;
	t20 = r_i_i_C(2) * t5;
	t18 = pkin(1) + r_i_i_C(3) + qJ(5) + pkin(8) + pkin(7);
	t16 = qJ(2) - t24;
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t4 = t14 * t20;
	t3 = t13 * t21;
	t2 = t22 + cos(qJ(3)) * pkin(3);
	t1 = [-t18 * t13 + t16 * t14, t13, t3 + (t2 - t20) * t13, t3 + (-t20 + t22) * t13, t14, 0; t16 * t13 + t18 * t14, -t14, t4 + (-t2 - t21) * t14, t4 + (-t21 - t22) * t14, t13, 0; 0, 0, t24, t15, 0, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:55
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (167->38), mult. (135->44), div. (0->0), fcn. (147->10), ass. (0->32)
	t21 = qJ(3) + qJ(4);
	t17 = pkin(10) + t21;
	t15 = sin(t17);
	t44 = pkin(9) + r_i_i_C(3);
	t46 = t44 * t15;
	t16 = cos(t17);
	t45 = t44 * t16 - pkin(4) * sin(t21);
	t42 = sin(qJ(3)) * pkin(3);
	t40 = pkin(4) * cos(t21);
	t22 = sin(qJ(6));
	t39 = r_i_i_C(2) * t22;
	t38 = pkin(1) + qJ(5) + pkin(8) + pkin(7);
	t26 = cos(qJ(1));
	t36 = t22 * t26;
	t24 = sin(qJ(1));
	t35 = t24 * t22;
	t25 = cos(qJ(6));
	t34 = t24 * t25;
	t33 = t25 * t26;
	t32 = t16 * t39;
	t31 = -r_i_i_C(1) * t25 - pkin(5);
	t30 = t24 * t46 + (pkin(5) * t24 + r_i_i_C(1) * t34) * t16;
	t29 = pkin(5) * t15 + qJ(2) + t42 - t45;
	t28 = t31 * t16 - t46;
	t27 = (t31 + t39) * t15 + t45;
	t8 = t40 + cos(qJ(3)) * pkin(3);
	t6 = t26 * t32;
	t4 = t15 * t33 - t35;
	t3 = t15 * t36 + t34;
	t2 = t15 * t34 + t36;
	t1 = -t15 * t35 + t33;
	t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t38 * t24 + t29 * t26, t24, (t8 - t32) * t24 + t30, (-t32 + t40) * t24 + t30, t26, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t29 * t24 + t38 * t26, -t26, t6 + (-t8 + t28) * t26, t6 + (t28 - t40) * t26, t24, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 0, t27 - t42, t27, 0, (-r_i_i_C(1) * t22 - r_i_i_C(2) * t25) * t16;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end