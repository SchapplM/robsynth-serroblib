% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRPPR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPPR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
Ja_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:52:14
	% EndTime: 2022-01-20 09:52:14
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:52:14
	% EndTime: 2022-01-20 09:52:14
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
	% StartTime: 2022-01-20 09:52:14
	% EndTime: 2022-01-20 09:52:14
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
	% StartTime: 2022-01-20 09:52:14
	% EndTime: 2022-01-20 09:52:14
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (30->8), mult. (14->8), div. (0->0), fcn. (14->6), ass. (0->7)
	t7 = qJ(1) + qJ(2);
	t5 = pkin(8) + t7;
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
	% StartTime: 2022-01-20 09:52:14
	% EndTime: 2022-01-20 09:52:15
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (71->12), mult. (34->10), div. (0->0), fcn. (36->8), ass. (0->9)
	t21 = qJ(4) + r_i_i_C(3);
	t20 = r_i_i_C(2) * sin(pkin(9)) - pkin(3) - r_i_i_C(1) * cos(pkin(9));
	t13 = qJ(1) + qJ(2);
	t11 = pkin(8) + t13;
	t8 = sin(t11);
	t9 = cos(t11);
	t17 = pkin(2) * cos(t13) + t21 * t8 - t20 * t9;
	t16 = -pkin(2) * sin(t13) + t21 * t9 + t20 * t8;
	t1 = [-sin(qJ(1)) * pkin(1) + t16, t16, 0, t8, 0; cos(qJ(1)) * pkin(1) + t17, t17, 0, -t9, 0; 0, 0, 1, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:52:14
	% EndTime: 2022-01-20 09:52:15
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (104->15), mult. (48->15), div. (0->0), fcn. (50->9), ass. (0->14)
	t24 = pkin(7) + qJ(4) + r_i_i_C(3);
	t14 = pkin(9) + qJ(5);
	t10 = sin(t14);
	t11 = cos(t14);
	t23 = t11 * r_i_i_C(1) - t10 * r_i_i_C(2);
	t22 = -cos(pkin(9)) * pkin(4) - pkin(3) - t23;
	t15 = qJ(1) + qJ(2);
	t19 = -r_i_i_C(1) * t10 - r_i_i_C(2) * t11;
	t12 = pkin(8) + t15;
	t6 = sin(t12);
	t7 = cos(t12);
	t18 = pkin(2) * cos(t15) + t24 * t6 - t22 * t7;
	t17 = -pkin(2) * sin(t15) + t24 * t7 + t22 * t6;
	t1 = [-sin(qJ(1)) * pkin(1) + t17, t17, 0, t6, t19 * t7; cos(qJ(1)) * pkin(1) + t18, t18, 0, -t7, t19 * t6; 0, 0, 1, 0, t23;];
	Ja_transl = t1;
end