% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:24
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPPR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
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
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
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
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (31->11), mult. (35->12), div. (0->0), fcn. (39->6), ass. (0->10)
	t12 = pkin(1) + r_i_i_C(3) + qJ(4) + pkin(7);
	t3 = qJ(3) + pkin(9);
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
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (67->18), mult. (76->18), div. (0->0), fcn. (87->8), ass. (0->14)
	t3 = qJ(3) + pkin(9);
	t1 = sin(t3);
	t4 = sin(pkin(10));
	t5 = cos(pkin(10));
	t13 = r_i_i_C(1) * t5 - r_i_i_C(2) * t4 + pkin(4);
	t14 = r_i_i_C(3) + qJ(5);
	t2 = cos(t3);
	t18 = -t13 * t1 + t14 * t2 - sin(qJ(3)) * pkin(3);
	t17 = t14 * t1 + t13 * t2 + cos(qJ(3)) * pkin(3);
	t12 = t4 * r_i_i_C(1) + t5 * r_i_i_C(2) + pkin(1) + pkin(7) + qJ(4);
	t11 = qJ(2) - t18;
	t10 = cos(qJ(1));
	t8 = sin(qJ(1));
	t6 = [t11 * t10 - t12 * t8, t8, t17 * t8, t10, -t8 * t2, 0; t12 * t10 + t11 * t8, -t10, -t17 * t10, t8, t10 * t2, 0; 0, 0, t18, 0, t1, 0;];
	Ja_transl = t6;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (116->31), mult. (105->38), div. (0->0), fcn. (120->10), ass. (0->24)
	t22 = r_i_i_C(3) + pkin(8) + qJ(5);
	t11 = qJ(3) + pkin(9);
	t9 = cos(t11);
	t30 = t22 * t9 - sin(qJ(3)) * pkin(3);
	t5 = cos(pkin(10)) * pkin(5) + pkin(4);
	t10 = pkin(10) + qJ(6);
	t6 = sin(t10);
	t8 = cos(t10);
	t20 = r_i_i_C(1) * t8 - r_i_i_C(2) * t6 + t5;
	t7 = sin(t11);
	t29 = t20 * t9 + t22 * t7 + cos(qJ(3)) * pkin(3);
	t16 = sin(qJ(1));
	t26 = t16 * t6;
	t25 = t16 * t8;
	t18 = cos(qJ(1));
	t24 = t18 * t6;
	t23 = t18 * t8;
	t21 = pkin(5) * sin(pkin(10)) + pkin(1) + qJ(4) + pkin(7);
	t19 = t7 * t5 + qJ(2) - t30;
	t4 = t7 * t23 - t26;
	t3 = t7 * t24 + t25;
	t2 = t7 * t25 + t24;
	t1 = -t7 * t26 + t23;
	t12 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t21 * t16 + t19 * t18, t16, t29 * t16, t18, -t16 * t9, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t19 * t16 + t21 * t18, -t18, -t29 * t18, t16, t18 * t9, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 0, -t20 * t7 + t30, 0, t7, (-r_i_i_C(1) * t6 - r_i_i_C(2) * t8) * t9;];
	Ja_transl = t12;
else
	Ja_transl=NaN(3,6);
end