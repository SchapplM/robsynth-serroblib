% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPRR4
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
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPPRR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPRR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
Ja_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
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
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (9->6), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t5 * t3 + t6 * t4, t3, 0, 0, 0; t6 * t3 + t5 * t4, -t4, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (18->11), mult. (30->12), div. (0->0), fcn. (39->6), ass. (0->9)
	t2 = sin(pkin(9));
	t3 = sin(pkin(8));
	t4 = cos(pkin(9));
	t5 = cos(pkin(8));
	t11 = pkin(1) - (-r_i_i_C(3) - qJ(3)) * t3 + (r_i_i_C(1) * t4 - r_i_i_C(2) * t2 + pkin(2)) * t5;
	t8 = t2 * r_i_i_C(1) + t4 * r_i_i_C(2) + qJ(2);
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [-t11 * t6 + t8 * t7, t6, t3 * t7, 0, 0; t11 * t7 + t8 * t6, -t7, t3 * t6, 0, 0; 0, 0, -t5, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (49->22), mult. (54->31), div. (0->0), fcn. (67->8), ass. (0->16)
	t11 = cos(pkin(8));
	t12 = sin(qJ(1));
	t16 = t11 * t12;
	t13 = cos(qJ(1));
	t15 = t11 * t13;
	t10 = sin(pkin(8));
	t14 = (cos(pkin(9)) * pkin(3) + pkin(2)) * t11 + pkin(1) + (r_i_i_C(3) + qJ(3) + pkin(6)) * t10;
	t9 = pkin(9) + qJ(4);
	t8 = cos(t9);
	t7 = sin(t9);
	t6 = sin(pkin(9)) * pkin(3) + qJ(2);
	t4 = t12 * t7 + t8 * t15;
	t3 = t12 * t8 - t7 * t15;
	t2 = t13 * t7 - t8 * t16;
	t1 = t13 * t8 + t7 * t16;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t14 * t12 + t6 * t13, t12, t10 * t13, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t6 * t12 + t14 * t13, -t13, t10 * t12, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0; 0, 0, -t11, (-r_i_i_C(1) * t7 - r_i_i_C(2) * t8) * t10, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (114->31), mult. (96->40), div. (0->0), fcn. (113->10), ass. (0->24)
	t17 = pkin(9) + qJ(4);
	t14 = cos(t17);
	t18 = sin(pkin(8));
	t19 = cos(pkin(8));
	t32 = (r_i_i_C(3) + pkin(7) + pkin(6) + qJ(3)) * t18 + t19 * (pkin(4) * t14 + cos(pkin(9)) * pkin(3) + pkin(2)) + pkin(1);
	t15 = qJ(5) + t17;
	t12 = cos(t15);
	t21 = cos(qJ(1));
	t11 = sin(t15);
	t20 = sin(qJ(1));
	t26 = t20 * t11;
	t5 = t12 * t21 + t19 * t26;
	t25 = t20 * t12;
	t6 = t11 * t21 - t19 * t25;
	t31 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t27 = t19 * t21;
	t7 = -t11 * t27 + t25;
	t8 = t12 * t27 + t26;
	t30 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t13 = sin(t17);
	t29 = pkin(4) * t13;
	t24 = qJ(2) + t29 + sin(pkin(9)) * pkin(3);
	t22 = -r_i_i_C(1) * t11 - r_i_i_C(2) * t12;
	t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t32 * t20 + t24 * t21, t20, t21 * t18, (-t13 * t27 + t14 * t20) * pkin(4) + t30, t30; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t24 * t20 + t32 * t21, -t21, t20 * t18, (-t13 * t19 * t20 - t14 * t21) * pkin(4) + t31, t31; 0, 0, -t19, (t22 - t29) * t18, t22 * t18;];
	Ja_transl = t1;
end