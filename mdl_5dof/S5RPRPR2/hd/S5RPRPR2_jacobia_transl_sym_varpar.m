% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR2
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
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRPR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:34:42
	% EndTime: 2020-01-03 11:34:42
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:34:42
	% EndTime: 2020-01-03 11:34:42
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:34:42
	% EndTime: 2020-01-03 11:34:42
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (8->5), mult. (6->6), div. (0->0), fcn. (6->4), ass. (0->4)
	t3 = qJ(1) + pkin(8);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [0, 1, 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + cos(qJ(1)) * pkin(1), 0, 0, 0, 0; t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + sin(qJ(1)) * pkin(1), 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:34:42
	% EndTime: 2020-01-03 11:34:42
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (26->8), mult. (12->8), div. (0->0), fcn. (12->6), ass. (0->7)
	t7 = qJ(1) + pkin(8);
	t6 = qJ(3) + t7;
	t4 = sin(t6);
	t5 = cos(t6);
	t9 = t4 * r_i_i_C(1) + t5 * r_i_i_C(2);
	t8 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t1 = [0, 1, 0, 0, 0; pkin(2) * cos(t7) + cos(qJ(1)) * pkin(1) + t8, 0, t8, 0, 0; pkin(2) * sin(t7) + sin(qJ(1)) * pkin(1) + t9, 0, t9, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:34:42
	% EndTime: 2020-01-03 11:34:42
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (68->13), mult. (32->10), div. (0->0), fcn. (34->8), ass. (0->9)
	t18 = r_i_i_C(3) + qJ(4);
	t17 = -r_i_i_C(2) * sin(pkin(9)) + r_i_i_C(1) * cos(pkin(9)) + pkin(3);
	t10 = qJ(1) + pkin(8);
	t9 = qJ(3) + t10;
	t7 = sin(t9);
	t8 = cos(t9);
	t14 = t17 * t8 + t18 * t7;
	t13 = t17 * t7 - t18 * t8;
	t1 = [0, 1, 0, 0, 0; pkin(2) * cos(t10) + cos(qJ(1)) * pkin(1) + t14, 0, t14, -t8, 0; pkin(2) * sin(t10) + sin(qJ(1)) * pkin(1) + t13, 0, t13, -t7, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:34:42
	% EndTime: 2020-01-03 11:34:42
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (101->17), mult. (46->15), div. (0->0), fcn. (48->9), ass. (0->14)
	t23 = r_i_i_C(3) + pkin(7) + qJ(4);
	t13 = pkin(9) + qJ(5);
	t10 = sin(t13);
	t11 = cos(t13);
	t22 = t11 * r_i_i_C(1) - t10 * r_i_i_C(2);
	t21 = cos(pkin(9)) * pkin(4) + pkin(3) + t22;
	t14 = qJ(1) + pkin(8);
	t18 = r_i_i_C(1) * t10 + r_i_i_C(2) * t11;
	t12 = qJ(3) + t14;
	t7 = sin(t12);
	t8 = cos(t12);
	t17 = t21 * t7 - t23 * t8;
	t16 = t21 * t8 + t23 * t7;
	t1 = [0, 1, 0, 0, t22; pkin(2) * cos(t14) + cos(qJ(1)) * pkin(1) + t16, 0, t16, -t8, -t18 * t7; pkin(2) * sin(t14) + sin(qJ(1)) * pkin(1) + t17, 0, t17, -t7, t18 * t8;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end