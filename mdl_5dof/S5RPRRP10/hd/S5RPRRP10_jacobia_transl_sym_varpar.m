% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRP10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRRP10_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRP10_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP10_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:52:35
	% EndTime: 2019-12-31 18:52:35
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:52:35
	% EndTime: 2019-12-31 18:52:35
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
	% StartTime: 2019-12-31 18:52:35
	% EndTime: 2019-12-31 18:52:35
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (9->6), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t3 * t5 + t4 * t6, t3, 0, 0, 0; t3 * t6 + t4 * t5, -t4, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:52:35
	% EndTime: 2019-12-31 18:52:35
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (26->9), mult. (26->11), div. (0->0), fcn. (28->5), ass. (0->10)
	t11 = r_i_i_C(3) + pkin(6) + qJ(2);
	t4 = pkin(8) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t10 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t9 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t8 = cos(pkin(8)) * pkin(2) + pkin(1) + t10;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [t11 * t7 - t8 * t6, t6, t9 * t7, 0, 0; t11 * t6 + t8 * t7, -t7, t9 * t6, 0, 0; 0, 0, t10, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:52:35
	% EndTime: 2019-12-31 18:52:35
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (71->24), mult. (85->33), div. (0->0), fcn. (95->7), ass. (0->22)
	t22 = pkin(7) + r_i_i_C(3);
	t8 = pkin(8) + qJ(3);
	t6 = sin(t8);
	t17 = t22 * t6;
	t7 = cos(t8);
	t23 = t17 + t7 * pkin(3) + cos(pkin(8)) * pkin(2) + pkin(1);
	t10 = sin(qJ(4));
	t11 = sin(qJ(1));
	t21 = t11 * t10;
	t12 = cos(qJ(4));
	t20 = t11 * t12;
	t13 = cos(qJ(1));
	t19 = t12 * t13;
	t18 = t13 * t10;
	t15 = r_i_i_C(1) * t12 - r_i_i_C(2) * t10 + pkin(3);
	t14 = -t15 * t6 + t22 * t7;
	t9 = -pkin(6) - qJ(2);
	t4 = t7 * t19 + t21;
	t3 = -t7 * t18 + t20;
	t2 = -t7 * t20 + t18;
	t1 = t7 * t21 + t19;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t23 * t11 - t13 * t9, t11, t14 * t13, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t11 * t9 + t23 * t13, -t13, t14 * t11, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, 0, t15 * t7 + t17, (-r_i_i_C(1) * t10 - r_i_i_C(2) * t12) * t6, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:52:35
	% EndTime: 2019-12-31 18:52:35
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (95->28), mult. (108->37), div. (0->0), fcn. (121->7), ass. (0->24)
	t25 = r_i_i_C(3) + qJ(5) + pkin(7);
	t9 = pkin(8) + qJ(3);
	t7 = sin(t9);
	t19 = t25 * t7;
	t14 = cos(qJ(4));
	t6 = pkin(4) * t14 + pkin(3);
	t8 = cos(t9);
	t27 = t19 + t8 * t6 + cos(pkin(8)) * pkin(2) + pkin(1);
	t26 = pkin(4) + r_i_i_C(1);
	t12 = sin(qJ(4));
	t15 = cos(qJ(1));
	t24 = t12 * t15;
	t13 = sin(qJ(1));
	t23 = t13 * t12;
	t22 = t13 * t14;
	t21 = t14 * t15;
	t18 = pkin(4) * t12 + pkin(6) + qJ(2);
	t17 = r_i_i_C(1) * t14 - r_i_i_C(2) * t12 + t6;
	t3 = -t8 * t24 + t22;
	t1 = t8 * t23 + t21;
	t16 = -t17 * t7 + t25 * t8;
	t4 = t8 * t21 + t23;
	t2 = -t8 * t22 + t24;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t27 * t13 + t18 * t15, t13, t16 * t15, -t4 * r_i_i_C(2) + t26 * t3, t15 * t7; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t18 * t13 + t27 * t15, -t15, t16 * t13, t2 * r_i_i_C(2) - t26 * t1, t13 * t7; 0, 0, t17 * t8 + t19, (-r_i_i_C(2) * t14 - t26 * t12) * t7, -t8;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,5);
end