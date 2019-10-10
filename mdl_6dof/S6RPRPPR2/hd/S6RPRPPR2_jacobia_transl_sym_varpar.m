% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:17
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPPR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
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
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (8->5), mult. (6->6), div. (0->0), fcn. (6->4), ass. (0->4)
	t3 = qJ(1) + pkin(9);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [-t1 * r_i_i_C(1) - t2 * r_i_i_C(2) - sin(qJ(1)) * pkin(1), 0, 0, 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + cos(qJ(1)) * pkin(1), 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (27->9), mult. (26->12), div. (0->0), fcn. (26->6), ass. (0->10)
	t9 = pkin(7) + r_i_i_C(3);
	t4 = sin(qJ(3));
	t5 = cos(qJ(3));
	t8 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t7 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5;
	t6 = pkin(2) + t8;
	t3 = qJ(1) + pkin(9);
	t2 = cos(t3);
	t1 = sin(t3);
	t10 = [-sin(qJ(1)) * pkin(1) + t9 * t2 - t6 * t1, 0, t7 * t2, 0, 0, 0; cos(qJ(1)) * pkin(1) + t9 * t1 + t6 * t2, 0, t7 * t1, 0, 0, 0; 0, 1, t8, 0, 0, 0;];
	Ja_transl = t10;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (49->13), mult. (33->14), div. (0->0), fcn. (35->8), ass. (0->11)
	t7 = qJ(3) + pkin(10);
	t2 = sin(t7);
	t4 = cos(t7);
	t15 = t4 * r_i_i_C(1) - t2 * r_i_i_C(2) + cos(qJ(3)) * pkin(3);
	t14 = r_i_i_C(3) + qJ(4) + pkin(7);
	t12 = pkin(2) + t15;
	t11 = -sin(qJ(3)) * pkin(3) - r_i_i_C(1) * t2 - r_i_i_C(2) * t4;
	t8 = qJ(1) + pkin(9);
	t5 = cos(t8);
	t3 = sin(t8);
	t1 = [-sin(qJ(1)) * pkin(1) + t14 * t5 - t12 * t3, 0, t11 * t5, t3, 0, 0; cos(qJ(1)) * pkin(1) + t14 * t3 + t12 * t5, 0, t11 * t3, -t5, 0, 0; 0, 1, t15, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (83->16), mult. (53->16), div. (0->0), fcn. (58->8), ass. (0->13)
	t13 = r_i_i_C(3) + qJ(5);
	t15 = pkin(4) - r_i_i_C(2);
	t7 = qJ(3) + pkin(10);
	t2 = sin(t7);
	t4 = cos(t7);
	t17 = cos(qJ(3)) * pkin(3) + t13 * t2 + t15 * t4;
	t16 = pkin(2) + t17;
	t14 = r_i_i_C(1) + qJ(4) + pkin(7);
	t11 = -sin(qJ(3)) * pkin(3) + t13 * t4 - t15 * t2;
	t8 = qJ(1) + pkin(9);
	t5 = cos(t8);
	t3 = sin(t8);
	t1 = [-sin(qJ(1)) * pkin(1) + t14 * t5 - t16 * t3, 0, t11 * t5, t3, t5 * t2, 0; cos(qJ(1)) * pkin(1) + t14 * t3 + t16 * t5, 0, t11 * t3, -t5, t3 * t2, 0; 0, 1, t17, 0, -t4, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (140->29), mult. (105->38), div. (0->0), fcn. (118->10), ass. (0->23)
	t21 = pkin(4) + pkin(8) + r_i_i_C(3);
	t11 = qJ(3) + pkin(10);
	t8 = cos(t11);
	t28 = cos(qJ(3)) * pkin(3) + t21 * t8;
	t6 = sin(t11);
	t27 = t6 * qJ(5) + pkin(2) + t28;
	t16 = cos(qJ(6));
	t12 = qJ(1) + pkin(9);
	t7 = sin(t12);
	t26 = t16 * t7;
	t9 = cos(t12);
	t25 = t16 * t9;
	t14 = sin(qJ(6));
	t24 = t7 * t14;
	t23 = t9 * t14;
	t22 = pkin(5) + qJ(4) + pkin(7);
	t18 = r_i_i_C(1) * t14 + r_i_i_C(2) * t16 + qJ(5);
	t17 = -sin(qJ(3)) * pkin(3) + t18 * t8 - t21 * t6;
	t4 = -t6 * t24 + t25;
	t3 = t6 * t26 + t23;
	t2 = t6 * t23 + t26;
	t1 = t6 * t25 - t24;
	t5 = [-sin(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t22 * t9 - t27 * t7, 0, t17 * t9, t7, t9 * t6, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; cos(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t22 * t7 + t27 * t9, 0, t17 * t7, -t9, t7 * t6, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 1, t18 * t6 + t28, 0, -t8, (-r_i_i_C(1) * t16 + r_i_i_C(2) * t14) * t8;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end