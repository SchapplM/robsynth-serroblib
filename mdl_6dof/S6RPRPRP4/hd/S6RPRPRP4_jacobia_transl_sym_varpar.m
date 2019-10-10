% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:34
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRP4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:14
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:14
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:14
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
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:14
	% DurationCPUTime: 0.13s
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
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:14
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (48->12), mult. (46->14), div. (0->0), fcn. (49->6), ass. (0->12)
	t10 = pkin(3) - r_i_i_C(2);
	t4 = sin(qJ(3));
	t5 = cos(qJ(3));
	t8 = r_i_i_C(3) + qJ(4);
	t7 = t10 * t5 + t8 * t4;
	t11 = pkin(2) + t7;
	t9 = pkin(7) + r_i_i_C(1);
	t6 = -t10 * t4 + t8 * t5;
	t3 = qJ(1) + pkin(9);
	t2 = cos(t3);
	t1 = sin(t3);
	t12 = [-sin(qJ(1)) * pkin(1) + t9 * t2 - t11 * t1, 0, t6 * t2, t2 * t4, 0, 0; cos(qJ(1)) * pkin(1) + t9 * t1 + t11 * t2, 0, t6 * t1, t1 * t4, 0, 0; 0, 1, t7, -t5, 0, 0;];
	Ja_transl = t12;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:14
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (89->25), mult. (98->38), div. (0->0), fcn. (109->8), ass. (0->20)
	t11 = cos(qJ(3));
	t16 = pkin(3) + pkin(8) + r_i_i_C(3);
	t14 = t16 * t11;
	t9 = sin(qJ(3));
	t20 = t9 * qJ(4) + pkin(2) + t14;
	t8 = sin(qJ(5));
	t19 = t8 * t9;
	t18 = pkin(4) + pkin(7);
	t10 = cos(qJ(5));
	t17 = t10 * t9;
	t13 = r_i_i_C(1) * t8 + r_i_i_C(2) * t10 + qJ(4);
	t12 = t13 * t11 - t16 * t9;
	t7 = qJ(1) + pkin(9);
	t6 = cos(t7);
	t5 = sin(t7);
	t4 = t10 * t6 - t5 * t19;
	t3 = t5 * t17 + t6 * t8;
	t2 = t10 * t5 + t6 * t19;
	t1 = t6 * t17 - t5 * t8;
	t15 = [-sin(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t18 * t6 - t20 * t5, 0, t12 * t6, t6 * t9, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; cos(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t18 * t5 + t20 * t6, 0, t12 * t5, t5 * t9, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0; 0, 1, t13 * t9 + t14, -t11, (-r_i_i_C(1) * t10 + r_i_i_C(2) * t8) * t11, 0;];
	Ja_transl = t15;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:14
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (139->28), mult. (157->38), div. (0->0), fcn. (181->8), ass. (0->23)
	t10 = sin(qJ(3));
	t12 = cos(qJ(3));
	t17 = pkin(3) + pkin(8) + r_i_i_C(2);
	t15 = t17 * t12;
	t24 = t10 * qJ(4) + pkin(2) + t15;
	t23 = pkin(4) + pkin(7);
	t22 = pkin(5) + r_i_i_C(1);
	t8 = qJ(1) + pkin(9);
	t6 = sin(t8);
	t21 = t6 * t10;
	t7 = cos(t8);
	t20 = t7 * t10;
	t11 = cos(qJ(5));
	t19 = t10 * t11;
	t18 = r_i_i_C(3) + qJ(6);
	t9 = sin(qJ(5));
	t14 = -t18 * t11 + t22 * t9 + qJ(4);
	t13 = -t17 * t10 + t14 * t12;
	t4 = t11 * t7 - t9 * t21;
	t3 = t6 * t19 + t7 * t9;
	t2 = t11 * t6 + t9 * t20;
	t1 = -t7 * t19 + t6 * t9;
	t5 = [-sin(qJ(1)) * pkin(1) + t23 * t7 + t22 * t4 + t18 * t3 - t24 * t6, 0, t13 * t7, t20, -t22 * t1 + t18 * t2, t1; cos(qJ(1)) * pkin(1) + t23 * t6 + t22 * t2 + t18 * t1 + t24 * t7, 0, t13 * t6, t21, -t18 * t4 + t22 * t3, -t3; 0, 1, t14 * t10 + t15, -t12, (-t22 * t11 - t18 * t9) * t12, t12 * t11;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end