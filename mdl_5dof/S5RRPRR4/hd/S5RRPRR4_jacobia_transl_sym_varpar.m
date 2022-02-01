% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR4
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
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRPRR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
Ja_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:49:11
	% EndTime: 2022-01-20 10:49:11
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:49:11
	% EndTime: 2022-01-20 10:49:11
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
	% StartTime: 2022-01-20 10:49:11
	% EndTime: 2022-01-20 10:49:11
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
	% StartTime: 2022-01-20 10:49:11
	% EndTime: 2022-01-20 10:49:11
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
	% StartTime: 2022-01-20 10:49:11
	% EndTime: 2022-01-20 10:49:11
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (77->12), mult. (44->14), div. (0->0), fcn. (44->8), ass. (0->13)
	t23 = r_i_i_C(3) + pkin(7);
	t14 = sin(qJ(4));
	t15 = cos(qJ(4));
	t22 = t15 * r_i_i_C(1) - t14 * r_i_i_C(2);
	t21 = -pkin(3) - t22;
	t13 = qJ(1) + qJ(2);
	t18 = -r_i_i_C(1) * t14 - r_i_i_C(2) * t15;
	t11 = pkin(9) + t13;
	t8 = sin(t11);
	t9 = cos(t11);
	t17 = pkin(2) * cos(t13) + t23 * t8 - t21 * t9;
	t16 = -pkin(2) * sin(t13) + t23 * t9 + t21 * t8;
	t1 = [-sin(qJ(1)) * pkin(1) + t16, t16, 0, t18 * t9, 0; cos(qJ(1)) * pkin(1) + t17, t17, 0, t18 * t8, 0; 0, 0, 1, t22, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:49:11
	% EndTime: 2022-01-20 10:49:11
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (123->15), mult. (63->18), div. (0->0), fcn. (63->10), ass. (0->16)
	t16 = qJ(4) + qJ(5);
	t12 = sin(t16);
	t14 = cos(t16);
	t24 = t14 * r_i_i_C(1) - r_i_i_C(2) * t12;
	t28 = cos(qJ(4)) * pkin(4) + t24;
	t27 = pkin(8) + pkin(7) + r_i_i_C(3);
	t26 = -pkin(3) - t28;
	t17 = qJ(1) + qJ(2);
	t23 = -r_i_i_C(1) * t12 - r_i_i_C(2) * t14;
	t22 = -sin(qJ(4)) * pkin(4) + t23;
	t11 = pkin(9) + t17;
	t6 = sin(t11);
	t7 = cos(t11);
	t21 = pkin(2) * cos(t17) + t27 * t6 - t26 * t7;
	t20 = -pkin(2) * sin(t17) + t27 * t7 + t26 * t6;
	t1 = [-sin(qJ(1)) * pkin(1) + t20, t20, 0, t22 * t7, t23 * t7; cos(qJ(1)) * pkin(1) + t21, t21, 0, t22 * t6, t23 * t6; 0, 0, 1, t28, t24;];
	Ja_transl = t1;
end