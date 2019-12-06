% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPPRP1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPRP1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:36:51
	% EndTime: 2019-12-05 17:36:51
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:36:51
	% EndTime: 2019-12-05 17:36:51
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
	% StartTime: 2019-12-05 17:36:51
	% EndTime: 2019-12-05 17:36:51
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (8->5), mult. (6->6), div. (0->0), fcn. (6->4), ass. (0->4)
	t3 = qJ(1) + pkin(7);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [0, 1, 0, 0, 0; -t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - cos(qJ(1)) * pkin(1), 0, 0, 0, 0; -t1 * r_i_i_C(1) - t2 * r_i_i_C(2) - sin(qJ(1)) * pkin(1), 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:36:51
	% EndTime: 2019-12-05 17:36:51
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (22->8), mult. (16->8), div. (0->0), fcn. (18->6), ass. (0->6)
	t7 = r_i_i_C(3) + qJ(3);
	t6 = -r_i_i_C(1) * cos(pkin(8)) + r_i_i_C(2) * sin(pkin(8)) - pkin(2);
	t3 = qJ(1) + pkin(7);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [0, 1, 0, 0, 0; -cos(qJ(1)) * pkin(1) - t7 * t1 + t6 * t2, 0, t2, 0, 0; -sin(qJ(1)) * pkin(1) + t7 * t2 + t6 * t1, 0, t1, 0, 0;];
	Ja_transl = t4;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:36:51
	% EndTime: 2019-12-05 17:36:51
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (53->19), mult. (54->27), div. (0->0), fcn. (64->8), ass. (0->17)
	t10 = sin(qJ(4));
	t7 = qJ(1) + pkin(7);
	t5 = sin(t7);
	t16 = t5 * t10;
	t11 = cos(qJ(4));
	t15 = t5 * t11;
	t6 = cos(t7);
	t14 = t6 * t10;
	t13 = t6 * t11;
	t8 = sin(pkin(8));
	t9 = cos(pkin(8));
	t12 = -pkin(3) * t9 - pkin(2) + (-pkin(6) - r_i_i_C(3)) * t8;
	t4 = -t9 * t13 - t16;
	t3 = t9 * t14 - t15;
	t2 = t9 * t15 - t14;
	t1 = t9 * t16 + t13;
	t17 = [0, 1, 0, (-r_i_i_C(1) * t10 - r_i_i_C(2) * t11) * t8, 0; -cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t5 * qJ(3) + t12 * t6, 0, t6, t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0; -sin(qJ(1)) * pkin(1) - t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t6 * qJ(3) + t12 * t5, 0, t5, -t3 * r_i_i_C(1) + t4 * r_i_i_C(2), 0;];
	Ja_transl = t17;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:36:51
	% EndTime: 2019-12-05 17:36:51
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (74->25), mult. (74->33), div. (0->0), fcn. (87->8), ass. (0->17)
	t18 = pkin(4) + r_i_i_C(1);
	t10 = cos(pkin(8));
	t12 = sin(qJ(4));
	t17 = t10 * t12;
	t13 = cos(qJ(4));
	t16 = t10 * t13;
	t15 = pkin(4) * t12 + qJ(3);
	t8 = qJ(1) + pkin(7);
	t6 = sin(t8);
	t7 = cos(t8);
	t3 = -t13 * t6 + t7 * t17;
	t1 = t13 * t7 + t6 * t17;
	t9 = sin(pkin(8));
	t14 = -t10 * (pkin(4) * t13 + pkin(3)) - pkin(2) + (-r_i_i_C(3) - qJ(5) - pkin(6)) * t9;
	t4 = -t6 * t12 - t7 * t16;
	t2 = -t7 * t12 + t6 * t16;
	t5 = [0, 1, 0, (-r_i_i_C(2) * t13 - t18 * t12) * t9, -t10; -cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t15 * t6 + t14 * t7, 0, t7, t2 * r_i_i_C(2) + t18 * t1, -t6 * t9; -sin(qJ(1)) * pkin(1) - t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * t7 + t14 * t6, 0, t6, t4 * r_i_i_C(2) - t18 * t3, t7 * t9;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,5);
end