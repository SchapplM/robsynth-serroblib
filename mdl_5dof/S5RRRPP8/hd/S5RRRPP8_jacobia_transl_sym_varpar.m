% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRPP8_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPP8_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP8_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_jacobia_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:10:08
	% EndTime: 2019-12-31 21:10:08
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:10:08
	% EndTime: 2019-12-31 21:10:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:10:08
	% EndTime: 2019-12-31 21:10:08
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (11->6), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->9)
	t8 = pkin(6) + r_i_i_C(3);
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t2 * t5 + t4 * t8, t6 * t4, 0, 0, 0; t2 * t8 + t4 * t5, t6 * t2, 0, 0, 0; 0, t7, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:10:08
	% EndTime: 2019-12-31 21:10:08
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (35->20), mult. (83->34), div. (0->0), fcn. (91->6), ass. (0->18)
	t16 = pkin(7) + r_i_i_C(3);
	t6 = sin(qJ(2));
	t14 = t16 * t6;
	t9 = cos(qJ(2));
	t18 = t9 * pkin(2) + pkin(1) + t14;
	t7 = sin(qJ(1));
	t17 = t7 * t9;
	t10 = cos(qJ(1));
	t15 = t10 * t9;
	t5 = sin(qJ(3));
	t8 = cos(qJ(3));
	t12 = r_i_i_C(1) * t8 - r_i_i_C(2) * t5 + pkin(2);
	t11 = -t12 * t6 + t16 * t9;
	t4 = t8 * t15 + t7 * t5;
	t3 = -t5 * t15 + t7 * t8;
	t2 = t10 * t5 - t8 * t17;
	t1 = t10 * t8 + t5 * t17;
	t13 = [t10 * pkin(6) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t18 * t7, t11 * t10, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0; t7 * pkin(6) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t18 * t10, t11 * t7, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0; 0, t12 * t9 + t14, (-r_i_i_C(1) * t5 - r_i_i_C(2) * t8) * t6, 0, 0;];
	Ja_transl = t13;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:10:08
	% EndTime: 2019-12-31 21:10:08
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (61->22), mult. (142->34), div. (0->0), fcn. (163->6), ass. (0->21)
	t11 = cos(qJ(2));
	t19 = pkin(7) + r_i_i_C(1);
	t8 = sin(qJ(2));
	t15 = t19 * t8;
	t23 = pkin(2) * t11 + pkin(1) + t15;
	t10 = cos(qJ(3));
	t16 = r_i_i_C(3) + qJ(4);
	t20 = pkin(3) - r_i_i_C(2);
	t7 = sin(qJ(3));
	t22 = t20 * t10 + t16 * t7 + pkin(2);
	t9 = sin(qJ(1));
	t21 = t9 * t7;
	t18 = t9 * t10;
	t12 = cos(qJ(1));
	t17 = t11 * t12;
	t13 = t19 * t11 - t22 * t8;
	t4 = t10 * t17 + t21;
	t3 = t7 * t17 - t18;
	t2 = t11 * t18 - t12 * t7;
	t1 = t10 * t12 + t11 * t21;
	t5 = [pkin(6) * t12 - t16 * t1 - t20 * t2 - t23 * t9, t13 * t12, t16 * t4 - t20 * t3, t3, 0; t9 * pkin(6) + t23 * t12 + t16 * t3 + t20 * t4, t13 * t9, -t20 * t1 + t16 * t2, t1, 0; 0, t22 * t11 + t15, (t16 * t10 - t20 * t7) * t8, t8 * t7, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:10:08
	% EndTime: 2019-12-31 21:10:08
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (80->22), mult. (184->35), div. (0->0), fcn. (214->6), ass. (0->21)
	t11 = cos(qJ(2));
	t17 = pkin(4) + pkin(7) + r_i_i_C(1);
	t8 = sin(qJ(2));
	t14 = t17 * t8;
	t23 = pkin(2) * t11 + pkin(1) + t14;
	t10 = cos(qJ(3));
	t16 = pkin(3) + r_i_i_C(3) + qJ(5);
	t18 = r_i_i_C(2) + qJ(4);
	t7 = sin(qJ(3));
	t22 = t16 * t10 + t18 * t7 + pkin(2);
	t9 = sin(qJ(1));
	t21 = t9 * t7;
	t20 = t9 * t10;
	t12 = cos(qJ(1));
	t19 = t11 * t12;
	t13 = t17 * t11 - t22 * t8;
	t4 = t10 * t19 + t21;
	t3 = t7 * t19 - t20;
	t2 = t11 * t20 - t12 * t7;
	t1 = t10 * t12 + t11 * t21;
	t5 = [pkin(6) * t12 - t18 * t1 - t16 * t2 - t23 * t9, t13 * t12, -t16 * t3 + t18 * t4, t3, t4; t9 * pkin(6) + t23 * t12 + t16 * t4 + t18 * t3, t13 * t9, -t16 * t1 + t18 * t2, t1, t2; 0, t22 * t11 + t14, (t18 * t10 - t16 * t7) * t8, t8 * t7, t8 * t10;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,5);
end