% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:01
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (23->9), mult. (16->8), div. (0->0), fcn. (18->6), ass. (0->6)
	t7 = r_i_i_C(3) + qJ(3);
	t6 = r_i_i_C(1) * cos(pkin(11)) - r_i_i_C(2) * sin(pkin(11)) + pkin(2);
	t3 = qJ(1) + pkin(10);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [-sin(qJ(1)) * pkin(1) + t7 * t2 - t6 * t1, 0, t1, 0, 0, 0; cos(qJ(1)) * pkin(1) + t7 * t1 + t6 * t2, 0, -t2, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (44->12), mult. (28->13), div. (0->0), fcn. (30->7), ass. (0->11)
	t12 = r_i_i_C(3) + pkin(7) + qJ(3);
	t6 = pkin(11) + qJ(4);
	t2 = sin(t6);
	t4 = cos(t6);
	t11 = r_i_i_C(1) * t4 - r_i_i_C(2) * t2;
	t10 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t4;
	t9 = cos(pkin(11)) * pkin(3) + pkin(2) + t11;
	t7 = qJ(1) + pkin(10);
	t5 = cos(t7);
	t3 = sin(t7);
	t1 = [-sin(qJ(1)) * pkin(1) + t12 * t5 - t9 * t3, 0, t3, t10 * t5, 0, 0; cos(qJ(1)) * pkin(1) + t12 * t3 + t9 * t5, 0, -t5, t10 * t3, 0, 0; 0, 1, 0, t11, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (87->15), mult. (45->17), div. (0->0), fcn. (47->9), ass. (0->14)
	t11 = pkin(11) + qJ(4);
	t9 = qJ(5) + t11;
	t4 = sin(t9);
	t5 = cos(t9);
	t16 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t19 = t16 + pkin(4) * cos(t11);
	t17 = r_i_i_C(3) + pkin(8) + pkin(7) + qJ(3);
	t15 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5;
	t14 = cos(pkin(11)) * pkin(3) + pkin(2) + t19;
	t13 = -pkin(4) * sin(t11) + t15;
	t12 = qJ(1) + pkin(10);
	t8 = cos(t12);
	t7 = sin(t12);
	t1 = [-sin(qJ(1)) * pkin(1) + t17 * t8 - t14 * t7, 0, t7, t13 * t8, t15 * t8, 0; cos(qJ(1)) * pkin(1) + t17 * t7 + t14 * t8, 0, -t8, t13 * t7, t15 * t7, 0; 0, 1, 0, t19, t16, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (208->35), mult. (125->42), div. (0->0), fcn. (135->11), ass. (0->31)
	t24 = pkin(11) + qJ(4);
	t22 = qJ(5) + t24;
	t17 = sin(t22);
	t18 = cos(t22);
	t26 = sin(qJ(6));
	t41 = r_i_i_C(2) * t26;
	t47 = pkin(9) + r_i_i_C(3);
	t48 = t17 * t41 + t18 * t47;
	t45 = t18 * pkin(5) + t47 * t17;
	t16 = pkin(4) * cos(t24);
	t44 = t16 + cos(pkin(11)) * pkin(3) + pkin(2) + t45;
	t27 = cos(qJ(6));
	t42 = r_i_i_C(1) * t27;
	t25 = qJ(1) + pkin(10);
	t20 = sin(t25);
	t38 = t20 * t26;
	t37 = t20 * t27;
	t21 = cos(t25);
	t36 = t21 * t26;
	t35 = t21 * t27;
	t34 = t48 * t20;
	t32 = t48 * t21;
	t30 = (-pkin(5) - t42) * t17;
	t29 = (-t41 + t42) * t18 + t45;
	t28 = -pkin(4) * sin(t24) + t30;
	t23 = -pkin(8) - pkin(7) - qJ(3);
	t4 = t18 * t35 + t38;
	t3 = -t18 * t36 + t37;
	t2 = -t18 * t37 + t36;
	t1 = t18 * t38 + t35;
	t5 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t21 * t23 - t44 * t20, 0, t20, t28 * t21 + t32, t21 * t30 + t32, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t20 * t23 + t44 * t21, 0, -t21, t28 * t20 + t34, t20 * t30 + t34, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, 1, 0, t16 + t29, t29, (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * t17;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end