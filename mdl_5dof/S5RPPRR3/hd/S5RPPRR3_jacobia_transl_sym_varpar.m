% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:40
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPPRR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPRR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:40:47
	% EndTime: 2019-10-24 10:40:47
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:40:47
	% EndTime: 2019-10-24 10:40:47
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
	% StartTime: 2019-10-24 10:40:47
	% EndTime: 2019-10-24 10:40:47
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
	% StartTime: 2019-10-24 10:40:47
	% EndTime: 2019-10-24 10:40:47
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (22->8), mult. (16->8), div. (0->0), fcn. (18->6), ass. (0->6)
	t7 = r_i_i_C(3) + qJ(3);
	t6 = -r_i_i_C(1) * cos(pkin(9)) + r_i_i_C(2) * sin(pkin(9)) - pkin(2);
	t3 = qJ(1) + pkin(8);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [0, 1, 0, 0, 0; -cos(qJ(1)) * pkin(1) - t7 * t1 + t6 * t2, 0, t2, 0, 0; -sin(qJ(1)) * pkin(1) + t7 * t2 + t6 * t1, 0, t1, 0, 0;];
	Ja_transl = t4;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:40:47
	% EndTime: 2019-10-24 10:40:47
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (43->12), mult. (28->13), div. (0->0), fcn. (30->7), ass. (0->11)
	t12 = r_i_i_C(3) + pkin(6) + qJ(3);
	t6 = pkin(9) + qJ(4);
	t2 = sin(t6);
	t4 = cos(t6);
	t11 = t4 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t10 = r_i_i_C(1) * t2 + r_i_i_C(2) * t4;
	t9 = -cos(pkin(9)) * pkin(3) - pkin(2) - t11;
	t7 = qJ(1) + pkin(8);
	t5 = cos(t7);
	t3 = sin(t7);
	t1 = [0, 1, 0, t11, 0; -cos(qJ(1)) * pkin(1) - t12 * t3 + t9 * t5, 0, t5, t10 * t3, 0; -sin(qJ(1)) * pkin(1) + t12 * t5 + t9 * t3, 0, t3, -t10 * t5, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:40:47
	% EndTime: 2019-10-24 10:40:47
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (86->16), mult. (45->17), div. (0->0), fcn. (47->9), ass. (0->15)
	t13 = pkin(9) + qJ(4);
	t11 = qJ(5) + t13;
	t6 = sin(t11);
	t7 = cos(t11);
	t25 = r_i_i_C(1) * t6 + r_i_i_C(2) * t7;
	t17 = t7 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t24 = t17 + pkin(4) * cos(t13);
	t23 = pkin(4) * sin(t13);
	t14 = qJ(1) + pkin(8);
	t9 = sin(t14);
	t19 = t25 * t9;
	t18 = r_i_i_C(3) + pkin(7) + pkin(6) + qJ(3);
	t15 = -cos(pkin(9)) * pkin(3) - pkin(2) - t24;
	t10 = cos(t14);
	t1 = [0, 1, 0, t24, t17; -cos(qJ(1)) * pkin(1) - t18 * t9 + t15 * t10, 0, t10, t9 * t23 + t19, t19; -sin(qJ(1)) * pkin(1) + t18 * t10 + t15 * t9, 0, t9, (-t25 - t23) * t10, -t25 * t10;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end