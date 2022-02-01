% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRPRR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
Ja_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:34:47
	% EndTime: 2022-01-20 10:34:48
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:34:48
	% EndTime: 2022-01-20 10:34:48
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
	% StartTime: 2022-01-20 10:34:47
	% EndTime: 2022-01-20 10:34:47
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (14->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->6)
	t4 = qJ(1) + qJ(2);
	t2 = sin(t4);
	t3 = cos(t4);
	t6 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t5 = -t2 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t1 = [-sin(qJ(1)) * pkin(1) + t5, t5, 0, 0, 0; cos(qJ(1)) * pkin(1) + t6, t6, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:34:47
	% EndTime: 2022-01-20 10:34:47
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (30->8), mult. (14->8), div. (0->0), fcn. (14->6), ass. (0->7)
	t7 = qJ(1) + qJ(2);
	t5 = pkin(9) + t7;
	t2 = sin(t5);
	t3 = cos(t5);
	t9 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2) + pkin(2) * cos(t7);
	t8 = -pkin(2) * sin(t7) - t2 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t1 = [-sin(qJ(1)) * pkin(1) + t8, t8, 0, 0, 0; cos(qJ(1)) * pkin(1) + t9, t9, 0, 0, 0; 0, 0, 1, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:34:47
	% EndTime: 2022-01-20 10:34:48
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (64->11), mult. (22->10), div. (0->0), fcn. (22->8), ass. (0->10)
	t10 = qJ(1) + qJ(2);
	t8 = pkin(9) + t10;
	t7 = qJ(4) + t8;
	t3 = sin(t7);
	t4 = cos(t7);
	t14 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t13 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t12 = t14 + pkin(3) * cos(t8) + pkin(2) * cos(t10);
	t11 = -pkin(2) * sin(t10) - pkin(3) * sin(t8) + t13;
	t1 = [-sin(qJ(1)) * pkin(1) + t11, t11, 0, t13, 0; cos(qJ(1)) * pkin(1) + t12, t12, 0, t14, 0; 0, 0, 1, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:34:47
	% EndTime: 2022-01-20 10:34:48
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (151->15), mult. (62->16), div. (0->0), fcn. (62->10), ass. (0->16)
	t28 = r_i_i_C(3) + pkin(8);
	t17 = sin(qJ(5));
	t18 = cos(qJ(5));
	t27 = t18 * r_i_i_C(1) - t17 * r_i_i_C(2);
	t26 = -pkin(4) - t27;
	t16 = qJ(1) + qJ(2);
	t14 = pkin(9) + t16;
	t23 = -r_i_i_C(1) * t17 - r_i_i_C(2) * t18;
	t13 = qJ(4) + t14;
	t10 = cos(t13);
	t9 = sin(t13);
	t22 = -t26 * t10 + t28 * t9;
	t21 = t28 * t10 + t26 * t9;
	t20 = pkin(2) * cos(t16) + t22 + pkin(3) * cos(t14);
	t19 = -pkin(2) * sin(t16) - pkin(3) * sin(t14) + t21;
	t1 = [-sin(qJ(1)) * pkin(1) + t19, t19, 0, t21, t23 * t10; cos(qJ(1)) * pkin(1) + t20, t20, 0, t22, t23 * t9; 0, 0, 1, 0, t27;];
	Ja_transl = t1;
end