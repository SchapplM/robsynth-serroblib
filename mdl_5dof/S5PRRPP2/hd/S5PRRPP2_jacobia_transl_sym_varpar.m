% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPP2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:29
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRPP2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPP2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPP2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:08
	% EndTime: 2019-10-24 10:29:09
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:08
	% EndTime: 2019-10-24 10:29:09
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:08
	% EndTime: 2019-10-24 10:29:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->2), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t3 = sin(qJ(2));
	t4 = cos(qJ(2));
	t5 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t1 = [0, t5 * cos(pkin(7)), 0, 0, 0; 0, t5 * sin(pkin(7)), 0, 0, 0; 1, r_i_i_C(1) * t4 - r_i_i_C(2) * t3, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:08
	% EndTime: 2019-10-24 10:29:09
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (19->12), mult. (51->25), div. (0->0), fcn. (55->6), ass. (0->12)
	t3 = sin(qJ(3));
	t6 = cos(qJ(2));
	t11 = t3 * t6;
	t5 = cos(qJ(3));
	t10 = t5 * t6;
	t9 = pkin(6) + r_i_i_C(3);
	t8 = t5 * r_i_i_C(1) - t3 * r_i_i_C(2) + pkin(2);
	t4 = sin(qJ(2));
	t7 = -t8 * t4 + t9 * t6;
	t2 = cos(pkin(7));
	t1 = sin(pkin(7));
	t12 = [0, t7 * t2, (t1 * t5 - t2 * t11) * r_i_i_C(1) + (-t1 * t3 - t2 * t10) * r_i_i_C(2), 0, 0; 0, t7 * t1, (-t1 * t11 - t2 * t5) * r_i_i_C(1) + (-t1 * t10 + t2 * t3) * r_i_i_C(2), 0, 0; 1, t9 * t4 + t8 * t6, (-r_i_i_C(1) * t3 - r_i_i_C(2) * t5) * t4, 0, 0;];
	Ja_transl = t12;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:08
	% EndTime: 2019-10-24 10:29:09
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (47->20), mult. (68->36), div. (0->0), fcn. (75->8), ass. (0->16)
	t17 = r_i_i_C(3) + qJ(4) + pkin(6);
	t11 = cos(qJ(2));
	t5 = sin(pkin(7));
	t16 = t11 * t5;
	t6 = cos(pkin(7));
	t15 = t11 * t6;
	t8 = sin(qJ(3));
	t14 = t11 * t8;
	t10 = cos(qJ(3));
	t4 = qJ(3) + pkin(8);
	t2 = sin(t4);
	t3 = cos(t4);
	t13 = pkin(3) * t10 + r_i_i_C(1) * t3 - r_i_i_C(2) * t2 + pkin(2);
	t9 = sin(qJ(2));
	t12 = t17 * t11 - t13 * t9;
	t1 = [0, t12 * t6, (-t2 * t15 + t5 * t3) * r_i_i_C(1) + (-t3 * t15 - t5 * t2) * r_i_i_C(2) + (t10 * t5 - t6 * t14) * pkin(3), t6 * t9, 0; 0, t12 * t5, (-t2 * t16 - t6 * t3) * r_i_i_C(1) + (-t3 * t16 + t6 * t2) * r_i_i_C(2) + (-t10 * t6 - t5 * t14) * pkin(3), t5 * t9, 0; 1, t13 * t11 + t17 * t9, (-pkin(3) * t8 - r_i_i_C(1) * t2 - r_i_i_C(2) * t3) * t9, -t11, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:08
	% EndTime: 2019-10-24 10:29:09
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (86->22), mult. (111->37), div. (0->0), fcn. (127->8), ass. (0->20)
	t15 = cos(qJ(3));
	t18 = r_i_i_C(3) + qJ(5);
	t23 = pkin(4) + r_i_i_C(1);
	t9 = qJ(3) + pkin(8);
	t7 = sin(t9);
	t8 = cos(t9);
	t24 = pkin(3) * t15 + t18 * t7 + t23 * t8 + pkin(2);
	t22 = r_i_i_C(2) + qJ(4) + pkin(6);
	t10 = sin(pkin(7));
	t16 = cos(qJ(2));
	t21 = t10 * t16;
	t11 = cos(pkin(7));
	t20 = t11 * t16;
	t13 = sin(qJ(3));
	t19 = t13 * t16;
	t14 = sin(qJ(2));
	t17 = -t24 * t14 + t22 * t16;
	t3 = -t10 * t8 + t7 * t20;
	t1 = t11 * t8 + t7 * t21;
	t2 = [0, t17 * t11, t18 * (t10 * t7 + t8 * t20) - t23 * t3 + (t10 * t15 - t11 * t19) * pkin(3), t11 * t14, t3; 0, t17 * t10, t18 * (-t11 * t7 + t8 * t21) - t23 * t1 + (-t10 * t19 - t11 * t15) * pkin(3), t10 * t14, t1; 1, t22 * t14 + t24 * t16, (-pkin(3) * t13 + t18 * t8 - t23 * t7) * t14, -t16, t14 * t7;];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,5);
end