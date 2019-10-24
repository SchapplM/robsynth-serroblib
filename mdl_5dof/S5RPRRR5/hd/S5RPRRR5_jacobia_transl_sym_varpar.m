% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRR5
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
% Datum: 2019-10-24 10:45
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRRR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:45:53
	% EndTime: 2019-10-24 10:45:53
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:45:53
	% EndTime: 2019-10-24 10:45:53
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
	% StartTime: 2019-10-24 10:45:53
	% EndTime: 2019-10-24 10:45:53
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
	% StartTime: 2019-10-24 10:45:53
	% EndTime: 2019-10-24 10:45:53
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (26->8), mult. (12->8), div. (0->0), fcn. (12->6), ass. (0->7)
	t5 = qJ(1) + pkin(9);
	t4 = qJ(3) + t5;
	t2 = sin(t4);
	t3 = cos(t4);
	t7 = -t3 * r_i_i_C(1) + t2 * r_i_i_C(2);
	t6 = -t2 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t1 = [0, 1, 0, 0, 0; -pkin(2) * cos(t5) - cos(qJ(1)) * pkin(1) + t7, 0, t7, 0, 0; -pkin(2) * sin(t5) - sin(qJ(1)) * pkin(1) + t6, 0, t6, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:45:53
	% EndTime: 2019-10-24 10:45:53
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (73->13), mult. (42->14), div. (0->0), fcn. (42->8), ass. (0->13)
	t10 = cos(qJ(4));
	t9 = sin(qJ(4));
	t19 = t10 * r_i_i_C(1) - t9 * r_i_i_C(2);
	t18 = pkin(7) + r_i_i_C(3);
	t17 = -pkin(3) - t19;
	t8 = qJ(1) + pkin(9);
	t13 = r_i_i_C(1) * t9 + r_i_i_C(2) * t10;
	t7 = qJ(3) + t8;
	t5 = sin(t7);
	t6 = cos(t7);
	t12 = t17 * t5 + t18 * t6;
	t11 = t17 * t6 - t18 * t5;
	t1 = [0, 1, 0, t19, 0; -pkin(2) * cos(t8) - cos(qJ(1)) * pkin(1) + t11, 0, t11, t13 * t5, 0; -pkin(2) * sin(t8) - sin(qJ(1)) * pkin(1) + t12, 0, t12, -t13 * t6, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:45:53
	% EndTime: 2019-10-24 10:45:53
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (119->19), mult. (61->21), div. (0->0), fcn. (61->10), ass. (0->21)
	t29 = r_i_i_C(3) + pkin(8) + pkin(7);
	t16 = qJ(4) + qJ(5);
	t12 = sin(t16);
	t15 = qJ(1) + pkin(9);
	t11 = qJ(3) + t15;
	t7 = sin(t11);
	t24 = t12 * t7;
	t13 = cos(t16);
	t26 = r_i_i_C(2) * t13;
	t28 = r_i_i_C(1) * t24 + t7 * t26;
	t27 = sin(qJ(4)) * pkin(4);
	t25 = t12 * r_i_i_C(2);
	t9 = t13 * r_i_i_C(1);
	t23 = t9 - t25;
	t14 = cos(qJ(4)) * pkin(4);
	t22 = -t14 - pkin(3) - t9;
	t21 = -r_i_i_C(1) * t12 - t26;
	t8 = cos(t11);
	t20 = (t22 + t25) * t8 - t29 * t7;
	t19 = r_i_i_C(2) * t24 + t22 * t7 + t29 * t8;
	t1 = [0, 1, 0, t14 + t23, t23; -pkin(2) * cos(t15) - cos(qJ(1)) * pkin(1) + t20, 0, t20, t7 * t27 + t28, t28; -pkin(2) * sin(t15) - sin(qJ(1)) * pkin(1) + t19, 0, t19, (t21 - t27) * t8, t21 * t8;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end