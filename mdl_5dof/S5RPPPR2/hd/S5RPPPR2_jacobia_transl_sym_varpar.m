% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPPPR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPPR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
Ja_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
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
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (9->6), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(7)) - r_i_i_C(2) * sin(pkin(7)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t5 * t3 + t6 * t4, t3, 0, 0, 0; t6 * t3 + t5 * t4, -t4, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (18->11), mult. (30->12), div. (0->0), fcn. (39->6), ass. (0->9)
	t2 = sin(pkin(8));
	t3 = sin(pkin(7));
	t4 = cos(pkin(8));
	t5 = cos(pkin(7));
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
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (34->24), mult. (61->35), div. (0->0), fcn. (81->8), ass. (0->21)
	t13 = sin(qJ(1));
	t8 = sin(pkin(8));
	t20 = t13 * t8;
	t14 = cos(qJ(1));
	t19 = t14 * t8;
	t9 = sin(pkin(7));
	t18 = t9 * t13;
	t17 = t9 * t14;
	t11 = cos(pkin(8));
	t16 = t13 * t11;
	t15 = t14 * t11;
	t12 = cos(pkin(7));
	t10 = cos(pkin(9));
	t7 = sin(pkin(9));
	t6 = -t8 * pkin(3) + qJ(4) * t11 - qJ(2);
	t5 = t12 * t15 + t20;
	t4 = t12 * t19 - t16;
	t3 = -t12 * t16 + t19;
	t2 = t12 * t20 + t15;
	t1 = (t11 * pkin(3) + t8 * qJ(4) + pkin(2)) * t12 + t9 * qJ(3) + pkin(1);
	t21 = [(t3 * t10 - t7 * t18) * r_i_i_C(1) + (-t10 * t18 - t3 * t7) * r_i_i_C(2) - t2 * r_i_i_C(3) - t1 * t13 - t6 * t14, t13, t17, t4, 0; (t5 * t10 + t7 * t17) * r_i_i_C(1) + (t10 * t17 - t5 * t7) * r_i_i_C(2) + t4 * r_i_i_C(3) + t1 * t14 - t6 * t13, -t14, t18, t2, 0; 0, 0, -t12, t9 * t8, 0;];
	Ja_transl = t21;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (79->45), mult. (163->78), div. (0->0), fcn. (213->10), ass. (0->33)
	t20 = sin(qJ(1));
	t22 = cos(qJ(1));
	t18 = cos(pkin(7));
	t21 = cos(qJ(5));
	t14 = sin(pkin(8));
	t19 = sin(qJ(5));
	t33 = t14 * t19;
	t13 = sin(pkin(9));
	t15 = sin(pkin(7));
	t16 = cos(pkin(9));
	t17 = cos(pkin(8));
	t29 = t16 * t17;
	t6 = t13 * t15 + t18 * t29;
	t24 = t18 * t33 + t21 * t6;
	t32 = t14 * t21;
	t9 = t16 * t32 - t17 * t19;
	t37 = -t24 * t20 + t22 * t9;
	t2 = t18 * t32 - t6 * t19;
	t7 = t16 * t33 + t17 * t21;
	t36 = t2 * t22 - t20 * t7;
	t31 = t15 * t14;
	t30 = t15 * t20;
	t28 = t16 * t22;
	t27 = t17 * t20;
	t26 = t18 * t22;
	t25 = t20 * t14;
	t11 = pkin(4) * t16 + pkin(6) * t13 + pkin(3);
	t23 = qJ(4) * t17 - t11 * t14 - qJ(2);
	t10 = t17 * t22 + t18 * t25;
	t8 = t14 * t26 - t27;
	t5 = -t13 * t18 + t15 * t29;
	t1 = (qJ(4) * t14 + t11 * t17 + pkin(2)) * t18 + pkin(1) + (t13 * pkin(4) - t16 * pkin(6) + qJ(3)) * t15;
	t3 = [t37 * r_i_i_C(1) + (-t2 * t20 - t22 * t7) * r_i_i_C(2) + ((t14 * t22 - t18 * t27) * t13 + t16 * t30) * r_i_i_C(3) - t1 * t20 - t23 * t22, t20, t15 * t22, t8, t36 * r_i_i_C(1) + (-t20 * t9 - t22 * t24) * r_i_i_C(2); ((t16 * t25 + t22 * t6) * t21 + t8 * t19) * r_i_i_C(1) + t36 * r_i_i_C(2) + ((t17 * t26 + t25) * t13 - t15 * t28) * r_i_i_C(3) + t1 * t22 - t23 * t20, -t22, t30, t10, (-(-t14 * t28 + t6 * t20) * t19 + t10 * t21) * r_i_i_C(1) + t37 * r_i_i_C(2); 0, 0, -t18, t31, (-t19 * t5 + t21 * t31) * r_i_i_C(1) + (-t19 * t31 - t21 * t5) * r_i_i_C(2);];
	Ja_transl = t3;
end