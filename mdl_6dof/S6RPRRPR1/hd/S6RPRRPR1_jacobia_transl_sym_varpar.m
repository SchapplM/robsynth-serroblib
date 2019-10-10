% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:23
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
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
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (8->5), mult. (6->6), div. (0->0), fcn. (6->4), ass. (0->4)
	t3 = qJ(1) + pkin(10);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [-t1 * r_i_i_C(1) - t2 * r_i_i_C(2) - sin(qJ(1)) * pkin(1), 0, 0, 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + cos(qJ(1)) * pkin(1), 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (27->9), mult. (26->12), div. (0->0), fcn. (26->6), ass. (0->10)
	t9 = pkin(7) + r_i_i_C(3);
	t4 = sin(qJ(3));
	t5 = cos(qJ(3));
	t8 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t7 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5;
	t6 = pkin(2) + t8;
	t3 = qJ(1) + pkin(10);
	t2 = cos(t3);
	t1 = sin(t3);
	t10 = [-sin(qJ(1)) * pkin(1) + t9 * t2 - t6 * t1, 0, t7 * t2, 0, 0, 0; cos(qJ(1)) * pkin(1) + t9 * t1 + t6 * t2, 0, t7 * t1, 0, 0, 0; 0, 1, t8, 0, 0, 0;];
	Ja_transl = t10;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (59->12), mult. (43->16), div. (0->0), fcn. (43->8), ass. (0->13)
	t9 = qJ(3) + qJ(4);
	t5 = sin(t9);
	t6 = cos(t9);
	t15 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t18 = t15 + cos(qJ(3)) * pkin(3);
	t16 = r_i_i_C(3) + pkin(8) + pkin(7);
	t14 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
	t13 = pkin(2) + t18;
	t12 = -sin(qJ(3)) * pkin(3) + t14;
	t8 = qJ(1) + pkin(10);
	t4 = cos(t8);
	t3 = sin(t8);
	t1 = [-sin(qJ(1)) * pkin(1) + t16 * t4 - t13 * t3, 0, t12 * t4, t14 * t4, 0, 0; cos(qJ(1)) * pkin(1) + t16 * t3 + t13 * t4, 0, t12 * t3, t14 * t3, 0, 0; 0, 1, t18, t15, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (98->16), mult. (53->18), div. (0->0), fcn. (55->10), ass. (0->14)
	t14 = qJ(3) + qJ(4);
	t9 = pkin(11) + t14;
	t4 = sin(t9);
	t5 = cos(t9);
	t19 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2) + pkin(4) * cos(t14);
	t23 = cos(qJ(3)) * pkin(3) + t19;
	t15 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5 - pkin(4) * sin(t14);
	t20 = r_i_i_C(3) + qJ(5) + pkin(8) + pkin(7);
	t17 = pkin(2) + t23;
	t16 = -sin(qJ(3)) * pkin(3) + t15;
	t13 = qJ(1) + pkin(10);
	t8 = cos(t13);
	t7 = sin(t13);
	t1 = [-sin(qJ(1)) * pkin(1) + t20 * t8 - t17 * t7, 0, t16 * t8, t15 * t8, t7, 0; cos(qJ(1)) * pkin(1) + t20 * t7 + t17 * t8, 0, t16 * t7, t15 * t7, -t8, 0; 0, 1, t23, t19, 0, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (219->35), mult. (133->43), div. (0->0), fcn. (143->12), ass. (0->31)
	t27 = qJ(3) + qJ(4);
	t22 = pkin(11) + t27;
	t17 = sin(t22);
	t18 = cos(t22);
	t28 = sin(qJ(6));
	t44 = r_i_i_C(2) * t28;
	t52 = pkin(9) + r_i_i_C(3);
	t53 = t17 * t44 + t18 * t52;
	t50 = t52 * t17 + t18 * pkin(5) + pkin(4) * cos(t27);
	t29 = cos(qJ(6));
	t45 = r_i_i_C(1) * t29;
	t31 = (-pkin(5) - t45) * t17 - pkin(4) * sin(t27);
	t24 = cos(qJ(3)) * pkin(3);
	t48 = pkin(2) + t24 + t50;
	t26 = qJ(1) + pkin(10);
	t20 = sin(t26);
	t41 = t20 * t28;
	t40 = t20 * t29;
	t21 = cos(t26);
	t39 = t21 * t28;
	t38 = t21 * t29;
	t37 = t53 * t20;
	t36 = t53 * t21;
	t32 = -sin(qJ(3)) * pkin(3) + t31;
	t30 = (-t44 + t45) * t18 + t50;
	t25 = -qJ(5) - pkin(8) - pkin(7);
	t4 = t18 * t38 + t41;
	t3 = -t18 * t39 + t40;
	t2 = -t18 * t40 + t39;
	t1 = t18 * t41 + t38;
	t5 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t21 * t25 - t48 * t20, 0, t32 * t21 + t36, t31 * t21 + t36, t20, t3 * r_i_i_C(1) - t4 * r_i_i_C(2); cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t20 * t25 + t48 * t21, 0, t32 * t20 + t37, t31 * t20 + t37, -t21, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2); 0, 1, t24 + t30, t30, 0, (-r_i_i_C(1) * t28 - r_i_i_C(2) * t29) * t17;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end