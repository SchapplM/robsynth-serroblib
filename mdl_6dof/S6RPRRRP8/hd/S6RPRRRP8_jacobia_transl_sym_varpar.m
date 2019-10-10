% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:55
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP8_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP8_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP8_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
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
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
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
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
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
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
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
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (103->33), mult. (125->42), div. (0->0), fcn. (135->8), ass. (0->30)
	t15 = qJ(3) + qJ(4);
	t13 = sin(t15);
	t39 = pkin(9) + r_i_i_C(3);
	t41 = t39 * t13;
	t14 = cos(t15);
	t40 = t39 * t14;
	t37 = sin(qJ(3)) * pkin(3);
	t36 = cos(qJ(3)) * pkin(3);
	t16 = sin(qJ(5));
	t35 = r_i_i_C(2) * t16;
	t34 = pkin(1) + pkin(8) + pkin(7);
	t18 = sin(qJ(1));
	t32 = t16 * t18;
	t21 = cos(qJ(1));
	t31 = t16 * t21;
	t19 = cos(qJ(5));
	t30 = t18 * t19;
	t29 = t19 * t21;
	t28 = t14 * t35;
	t27 = -r_i_i_C(1) * t19 - pkin(4);
	t26 = t18 * t41 + (pkin(4) * t18 + r_i_i_C(1) * t30) * t14;
	t25 = t40 + (t27 + t35) * t13;
	t24 = t27 * t14 - t41;
	t23 = pkin(4) * t13 + qJ(2) + t37 - t40;
	t6 = t21 * t28;
	t4 = t13 * t29 - t32;
	t3 = t13 * t31 + t30;
	t2 = t13 * t30 + t31;
	t1 = -t13 * t32 + t29;
	t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t34 * t18 + t23 * t21, t18, (-t28 + t36) * t18 + t26, -t18 * t28 + t26, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t23 * t18 + t34 * t21, -t21, t6 + (t24 - t36) * t21, t24 * t21 + t6, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0; 0, 0, t25 - t37, t25, (-r_i_i_C(1) * t16 - r_i_i_C(2) * t19) * t14, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (160->33), mult. (200->41), div. (0->0), fcn. (223->8), ass. (0->29)
	t17 = qJ(3) + qJ(4);
	t15 = sin(t17);
	t43 = pkin(9) + r_i_i_C(2);
	t45 = t15 * t43;
	t16 = cos(t17);
	t44 = t43 * t16;
	t42 = pkin(5) + r_i_i_C(1);
	t32 = r_i_i_C(3) + qJ(6);
	t40 = cos(qJ(3)) * pkin(3);
	t39 = sin(qJ(3)) * pkin(3);
	t38 = pkin(1) + pkin(8) + pkin(7);
	t18 = sin(qJ(5));
	t20 = sin(qJ(1));
	t36 = t20 * t18;
	t21 = cos(qJ(5));
	t35 = t20 * t21;
	t23 = cos(qJ(1));
	t34 = t23 * t18;
	t33 = t23 * t21;
	t29 = t20 * t45 + (pkin(4) * t20 + t32 * t36 + t35 * t42) * t16;
	t28 = -t32 * t18 - t42 * t21 - pkin(4);
	t27 = t15 * pkin(4) + qJ(2) + t39 - t44;
	t26 = t28 * t15 + t44;
	t25 = t28 * t16 - t45;
	t4 = t15 * t33 - t36;
	t3 = t15 * t34 + t35;
	t2 = t15 * t35 + t34;
	t1 = t15 * t36 - t33;
	t5 = [-t38 * t20 + t27 * t23 + t32 * t3 + t42 * t4, t20, t20 * t40 + t29, t29, -t42 * t1 + t32 * t2, t1; t32 * t1 + t42 * t2 + t27 * t20 + t38 * t23, -t23, (t25 - t40) * t23, t25 * t23, t42 * t3 - t32 * t4, -t3; 0, 0, t26 - t39, t26, (-t42 * t18 + t32 * t21) * t16, t16 * t18;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end