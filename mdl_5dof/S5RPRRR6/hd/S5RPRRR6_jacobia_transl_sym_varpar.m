% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:46
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRRR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:46:15
	% EndTime: 2019-10-24 10:46:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:46:15
	% EndTime: 2019-10-24 10:46:15
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
	% StartTime: 2019-10-24 10:46:15
	% EndTime: 2019-10-24 10:46:15
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (8->5), mult. (6->6), div. (0->0), fcn. (6->4), ass. (0->4)
	t3 = qJ(1) + pkin(9);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [0, 1, 0, 0, 0; -t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - cos(qJ(1)) * pkin(1), 0, 0, 0, 0; -t1 * r_i_i_C(1) - t2 * r_i_i_C(2) - sin(qJ(1)) * pkin(1), 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:46:15
	% EndTime: 2019-10-24 10:46:15
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (27->10), mult. (26->12), div. (0->0), fcn. (26->6), ass. (0->10)
	t9 = pkin(6) + r_i_i_C(3);
	t4 = sin(qJ(3));
	t5 = cos(qJ(3));
	t8 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t7 = r_i_i_C(1) * t4 + r_i_i_C(2) * t5;
	t6 = -pkin(2) - t8;
	t3 = qJ(1) + pkin(9);
	t2 = cos(t3);
	t1 = sin(t3);
	t10 = [0, 1, t8, 0, 0; -cos(qJ(1)) * pkin(1) - t9 * t1 + t6 * t2, 0, t7 * t1, 0, 0; -sin(qJ(1)) * pkin(1) + t9 * t2 + t6 * t1, 0, -t7 * t2, 0, 0;];
	Ja_transl = t10;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:46:15
	% EndTime: 2019-10-24 10:46:15
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (59->14), mult. (43->16), div. (0->0), fcn. (43->8), ass. (0->14)
	t11 = qJ(3) + qJ(4);
	t7 = sin(t11);
	t8 = cos(t11);
	t24 = r_i_i_C(1) * t7 + r_i_i_C(2) * t8;
	t16 = t8 * r_i_i_C(1) - t7 * r_i_i_C(2);
	t23 = t16 + cos(qJ(3)) * pkin(3);
	t10 = qJ(1) + pkin(9);
	t5 = sin(t10);
	t19 = t24 * t5;
	t18 = sin(qJ(3)) * pkin(3);
	t17 = r_i_i_C(3) + pkin(7) + pkin(6);
	t14 = -pkin(2) - t23;
	t6 = cos(t10);
	t1 = [0, 1, t23, t16, 0; -cos(qJ(1)) * pkin(1) - t17 * t5 + t14 * t6, 0, t18 * t5 + t19, t19, 0; -sin(qJ(1)) * pkin(1) + t17 * t6 + t14 * t5, 0, (-t24 - t18) * t6, -t24 * t6, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:46:15
	% EndTime: 2019-10-24 10:46:15
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (114->19), mult. (63->20), div. (0->0), fcn. (63->10), ass. (0->17)
	t16 = qJ(3) + qJ(4);
	t12 = qJ(5) + t16;
	t7 = sin(t12);
	t8 = cos(t12);
	t20 = t8 * r_i_i_C(1) - t7 * r_i_i_C(2);
	t19 = t20 + pkin(4) * cos(t16);
	t28 = cos(qJ(3)) * pkin(3) + t19;
	t27 = r_i_i_C(1) * t7 + r_i_i_C(2) * t8;
	t14 = qJ(1) + pkin(9);
	t9 = sin(t14);
	t23 = t27 * t9;
	t22 = pkin(4) * sin(t16);
	t21 = r_i_i_C(3) + pkin(8) + pkin(7) + pkin(6);
	t17 = -pkin(2) - t28;
	t10 = cos(t14);
	t4 = -t22 - sin(qJ(3)) * pkin(3);
	t1 = [0, 1, t28, t19, t20; -cos(qJ(1)) * pkin(1) - t21 * t9 + t17 * t10, 0, -t9 * t4 + t23, t9 * t22 + t23, t23; -sin(qJ(1)) * pkin(1) + t21 * t10 + t17 * t9, 0, (-t27 + t4) * t10, (-t27 - t22) * t10, -t27 * t10;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end