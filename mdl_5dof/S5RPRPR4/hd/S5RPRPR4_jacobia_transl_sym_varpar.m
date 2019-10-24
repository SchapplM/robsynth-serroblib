% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:42
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRPR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:35
	% EndTime: 2019-10-24 10:42:35
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:35
	% EndTime: 2019-10-24 10:42:35
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [0, 0, 0, 0, 0; -r_i_i_C(1) * t2 + r_i_i_C(2) * t1, 0, 0, 0, 0; -r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:35
	% EndTime: 2019-10-24 10:42:35
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (8->5), mult. (6->6), div. (0->0), fcn. (6->4), ass. (0->4)
	t3 = qJ(1) + pkin(8);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [0, 1, 0, 0, 0; -t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - cos(qJ(1)) * pkin(1), 0, 0, 0, 0; -t1 * r_i_i_C(1) - t2 * r_i_i_C(2) - sin(qJ(1)) * pkin(1), 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:35
	% EndTime: 2019-10-24 10:42:35
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (27->10), mult. (26->12), div. (0->0), fcn. (26->6), ass. (0->10)
	t9 = pkin(6) + r_i_i_C(3);
	t4 = sin(qJ(3));
	t5 = cos(qJ(3));
	t8 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t7 = r_i_i_C(1) * t4 + r_i_i_C(2) * t5;
	t6 = -pkin(2) - t8;
	t3 = qJ(1) + pkin(8);
	t2 = cos(t3);
	t1 = sin(t3);
	t10 = [0, 1, t8, 0, 0; -cos(qJ(1)) * pkin(1) - t9 * t1 + t6 * t2, 0, t7 * t1, 0, 0; -sin(qJ(1)) * pkin(1) + t9 * t2 + t6 * t1, 0, -t7 * t2, 0, 0;];
	Ja_transl = t10;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:35
	% EndTime: 2019-10-24 10:42:35
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (48->13), mult. (33->14), div. (0->0), fcn. (35->8), ass. (0->11)
	t7 = qJ(3) + pkin(9);
	t2 = sin(t7);
	t4 = cos(t7);
	t15 = t4 * r_i_i_C(1) - t2 * r_i_i_C(2) + cos(qJ(3)) * pkin(3);
	t14 = r_i_i_C(3) + qJ(4) + pkin(6);
	t12 = -pkin(2) - t15;
	t11 = sin(qJ(3)) * pkin(3) + r_i_i_C(1) * t2 + r_i_i_C(2) * t4;
	t8 = qJ(1) + pkin(8);
	t5 = cos(t8);
	t3 = sin(t8);
	t1 = [0, 1, t15, 0, 0; -cos(qJ(1)) * pkin(1) - t14 * t3 + t12 * t5, 0, t11 * t3, t5, 0; -sin(qJ(1)) * pkin(1) + t14 * t5 + t12 * t3, 0, -t11 * t5, t3, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:35
	% EndTime: 2019-10-24 10:42:35
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (89->17), mult. (48->18), div. (0->0), fcn. (50->10), ass. (0->15)
	t14 = qJ(3) + pkin(9);
	t11 = qJ(5) + t14;
	t7 = sin(t11);
	t8 = cos(t11);
	t18 = t8 * r_i_i_C(1) - t7 * r_i_i_C(2);
	t26 = t18 + cos(qJ(3)) * pkin(3) + pkin(4) * cos(t14);
	t25 = r_i_i_C(1) * t7 + r_i_i_C(2) * t8;
	t15 = qJ(1) + pkin(8);
	t9 = sin(t15);
	t21 = t25 * t9;
	t20 = r_i_i_C(3) + pkin(7) + qJ(4) + pkin(6);
	t16 = -pkin(2) - t26;
	t10 = cos(t15);
	t4 = -pkin(4) * sin(t14) - sin(qJ(3)) * pkin(3);
	t1 = [0, 1, t26, 0, t18; -cos(qJ(1)) * pkin(1) - t20 * t9 + t16 * t10, 0, -t9 * t4 + t21, t10, t21; -sin(qJ(1)) * pkin(1) + t20 * t10 + t16 * t9, 0, (-t25 + t4) * t10, t9, -t25 * t10;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end