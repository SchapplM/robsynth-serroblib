% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRRP2
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:19
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PPRRP2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRRP2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:40
	% EndTime: 2019-10-24 10:19:40
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:40
	% EndTime: 2019-10-24 10:19:41
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:40
	% EndTime: 2019-10-24 10:19:40
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(pkin(7)), 0, 0, 0; 0, -cos(pkin(7)), 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:40
	% EndTime: 2019-10-24 10:19:40
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (10->4), mult. (10->6), div. (0->0), fcn. (12->4), ass. (0->7)
	t3 = pkin(8) + qJ(3);
	t1 = sin(t3);
	t2 = cos(t3);
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t2;
	t5 = cos(pkin(7));
	t4 = sin(pkin(7));
	t7 = [0, t4, t6 * t5, 0, 0; 0, -t5, t6 * t4, 0, 0; 1, 0, t2 * r_i_i_C(1) - t1 * r_i_i_C(2), 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:40
	% EndTime: 2019-10-24 10:19:40
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (41->14), mult. (51->23), div. (0->0), fcn. (57->6), ass. (0->15)
	t4 = sin(pkin(7));
	t6 = sin(qJ(4));
	t14 = t4 * t6;
	t7 = cos(qJ(4));
	t13 = t4 * t7;
	t5 = cos(pkin(7));
	t12 = t5 * t6;
	t11 = t5 * t7;
	t10 = pkin(6) + r_i_i_C(3);
	t9 = t7 * r_i_i_C(1) - t6 * r_i_i_C(2) + pkin(3);
	t3 = pkin(8) + qJ(3);
	t1 = sin(t3);
	t2 = cos(t3);
	t8 = -t9 * t1 + t10 * t2;
	t15 = [0, t4, t8 * t5, (-t2 * t12 + t13) * r_i_i_C(1) + (-t2 * t11 - t14) * r_i_i_C(2), 0; 0, -t5, t8 * t4, (-t2 * t14 - t11) * r_i_i_C(1) + (-t2 * t13 + t12) * r_i_i_C(2), 0; 1, 0, t10 * t1 + t9 * t2, (-r_i_i_C(1) * t6 - r_i_i_C(2) * t7) * t1, 0;];
	Ja_transl = t15;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:40
	% EndTime: 2019-10-24 10:19:41
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (74->16), mult. (94->24), div. (0->0), fcn. (109->6), ass. (0->19)
	t11 = sin(qJ(4));
	t12 = cos(qJ(4));
	t14 = r_i_i_C(3) + qJ(5);
	t20 = pkin(4) + r_i_i_C(1);
	t21 = t14 * t11 + t20 * t12 + pkin(3);
	t19 = pkin(6) + r_i_i_C(2);
	t9 = sin(pkin(7));
	t18 = t9 * t11;
	t17 = t9 * t12;
	t10 = cos(pkin(7));
	t16 = t10 * t11;
	t15 = t10 * t12;
	t8 = pkin(8) + qJ(3);
	t6 = sin(t8);
	t7 = cos(t8);
	t13 = t19 * t7 - t21 * t6;
	t3 = t7 * t16 - t17;
	t1 = t7 * t18 + t15;
	t2 = [0, t9, t13 * t10, t14 * (t7 * t15 + t18) - t20 * t3, t3; 0, -t10, t13 * t9, t14 * (t7 * t17 - t16) - t20 * t1, t1; 1, 0, t19 * t6 + t21 * t7, (-t20 * t11 + t14 * t12) * t6, t6 * t11;];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,5);
end