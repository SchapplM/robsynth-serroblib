% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRPR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:10:13
	% EndTime: 2020-01-03 12:10:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:10:13
	% EndTime: 2020-01-03 12:10:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:10:13
	% EndTime: 2020-01-03 12:10:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (14->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->6)
	t6 = qJ(1) + qJ(2);
	t4 = sin(t6);
	t5 = cos(t6);
	t8 = t4 * r_i_i_C(1) + t5 * r_i_i_C(2);
	t7 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t1 = [0, 0, 0, 0, 0; cos(qJ(1)) * pkin(1) + t7, t7, 0, 0, 0; sin(qJ(1)) * pkin(1) + t8, t8, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:10:13
	% EndTime: 2020-01-03 12:10:13
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (45->10), mult. (40->12), div. (0->0), fcn. (40->6), ass. (0->12)
	t19 = pkin(7) + r_i_i_C(3);
	t10 = sin(qJ(3));
	t11 = cos(qJ(3));
	t18 = t11 * r_i_i_C(1) - t10 * r_i_i_C(2);
	t17 = pkin(2) + t18;
	t14 = r_i_i_C(1) * t10 + r_i_i_C(2) * t11;
	t9 = qJ(1) + qJ(2);
	t7 = sin(t9);
	t8 = cos(t9);
	t13 = t17 * t8 + t19 * t7;
	t12 = t17 * t7 - t19 * t8;
	t1 = [0, 0, t18, 0, 0; cos(qJ(1)) * pkin(1) + t13, t13, -t14 * t7, 0, 0; sin(qJ(1)) * pkin(1) + t12, t12, t14 * t8, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:10:13
	% EndTime: 2020-01-03 12:10:13
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (76->15), mult. (49->14), div. (0->0), fcn. (51->8), ass. (0->13)
	t13 = qJ(3) + pkin(9);
	t8 = sin(t13);
	t9 = cos(t13);
	t25 = cos(qJ(3)) * pkin(3) + t9 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t24 = r_i_i_C(3) + qJ(4) + pkin(7);
	t22 = pkin(2) + t25;
	t19 = sin(qJ(3)) * pkin(3) + r_i_i_C(1) * t8 + r_i_i_C(2) * t9;
	t14 = qJ(1) + qJ(2);
	t10 = sin(t14);
	t11 = cos(t14);
	t18 = t22 * t10 - t24 * t11;
	t17 = t24 * t10 + t22 * t11;
	t1 = [0, 0, t25, 0, 0; cos(qJ(1)) * pkin(1) + t17, t17, -t19 * t10, -t11, 0; sin(qJ(1)) * pkin(1) + t18, t18, t19 * t11, -t10, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:10:13
	% EndTime: 2020-01-03 12:10:13
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (127->22), mult. (66->23), div. (0->0), fcn. (68->10), ass. (0->19)
	t31 = r_i_i_C(3) + pkin(8) + qJ(4) + pkin(7);
	t20 = qJ(3) + pkin(9);
	t15 = qJ(5) + t20;
	t12 = sin(t15);
	t13 = cos(t15);
	t25 = t13 * r_i_i_C(1) - t12 * r_i_i_C(2);
	t21 = qJ(1) + qJ(2);
	t17 = cos(t21);
	t27 = t13 * t17;
	t28 = t12 * t17;
	t30 = r_i_i_C(1) * t28 + r_i_i_C(2) * t27;
	t26 = pkin(4) * cos(t20) + cos(qJ(3)) * pkin(3);
	t24 = -r_i_i_C(1) * t12 - r_i_i_C(2) * t13;
	t16 = sin(t21);
	t3 = pkin(2) + t26;
	t23 = -t31 * t17 + (t3 + t25) * t16;
	t22 = r_i_i_C(1) * t27 - r_i_i_C(2) * t28 + t31 * t16 + t17 * t3;
	t8 = -pkin(4) * sin(t20) - sin(qJ(3)) * pkin(3);
	t1 = [0, 0, t25 + t26, 0, t25; cos(qJ(1)) * pkin(1) + t22, t22, (t24 + t8) * t16, -t17, t24 * t16; sin(qJ(1)) * pkin(1) + t23, t23, -t17 * t8 + t30, -t16, t30;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end