% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRR13
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRRR13_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR13_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR13_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
Ja_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (14->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->6)
	t4 = qJ(1) + qJ(2);
	t2 = sin(t4);
	t3 = cos(t4);
	t6 = r_i_i_C(1) * t3 - r_i_i_C(2) * t2;
	t5 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t1 = [-sin(qJ(1)) * pkin(1) + t5, t5, 0, 0, 0; cos(qJ(1)) * pkin(1) + t6, t6, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (40->8), mult. (18->8), div. (0->0), fcn. (18->6), ass. (0->9)
	t7 = qJ(1) + qJ(2);
	t6 = qJ(3) + t7;
	t3 = sin(t6);
	t4 = cos(t6);
	t11 = r_i_i_C(1) * t4 - r_i_i_C(2) * t3;
	t10 = t11 + pkin(2) * cos(t7);
	t9 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t8 = -pkin(2) * sin(t7) + t9;
	t1 = [-sin(qJ(1)) * pkin(1) + t8, t8, t9, 0, 0; cos(qJ(1)) * pkin(1) + t10, t10, t11, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (153->20), mult. (104->30), div. (0->0), fcn. (120->10), ass. (0->20)
	t20 = sin(pkin(5));
	t32 = t20 * (pkin(9) + r_i_i_C(3));
	t21 = cos(pkin(5));
	t22 = sin(qJ(4));
	t29 = t21 * t22;
	t23 = cos(qJ(4));
	t28 = t21 * t23;
	t19 = qJ(1) + qJ(2);
	t18 = qJ(3) + t19;
	t15 = sin(t18);
	t16 = cos(t18);
	t7 = -t15 * t28 - t16 * t22;
	t8 = -t15 * t29 + t16 * t23;
	t27 = t16 * pkin(3) + t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t15 * t32;
	t26 = pkin(2) * cos(t19) + t27;
	t5 = t15 * t22 - t16 * t28;
	t6 = -t15 * t23 - t16 * t29;
	t25 = -t15 * pkin(3) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t16 * t32;
	t24 = -pkin(2) * sin(t19) + t25;
	t1 = [-sin(qJ(1)) * pkin(1) + t24, t24, t25, t7 * r_i_i_C(1) - t8 * r_i_i_C(2), 0; cos(qJ(1)) * pkin(1) + t26, t26, t27, -t5 * r_i_i_C(1) + t6 * r_i_i_C(2), 0; 0, 0, 0, (r_i_i_C(1) * t23 - r_i_i_C(2) * t22) * t20, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (347->35), mult. (178->44), div. (0->0), fcn. (174->16), ass. (0->35)
	t34 = sin(pkin(5));
	t35 = cos(pkin(5));
	t36 = sin(qJ(4));
	t53 = -t35 * t36 * pkin(4) + (pkin(9) + pkin(10) + r_i_i_C(3)) * t34;
	t32 = qJ(4) + qJ(5);
	t45 = pkin(5) - t32;
	t44 = sin(t45);
	t27 = pkin(5) + t32;
	t46 = sin(t27) / 0.2e1;
	t13 = t46 - t44 / 0.2e1;
	t33 = qJ(1) + qJ(2);
	t31 = qJ(3) + t33;
	t24 = sin(t31);
	t25 = cos(t31);
	t30 = cos(t32);
	t43 = -t25 * t13 - t24 * t30;
	t21 = cos(t27) / 0.2e1;
	t22 = cos(t45);
	t14 = t22 / 0.2e1 + t21;
	t28 = sin(t32);
	t9 = -t25 * t14 + t24 * t28;
	t52 = -t9 * r_i_i_C(1) + t43 * r_i_i_C(2);
	t10 = -t24 * t14 - t25 * t28;
	t42 = t24 * t13 - t25 * t30;
	t51 = t10 * r_i_i_C(1) + t42 * r_i_i_C(2);
	t37 = cos(qJ(4));
	t49 = t37 * pkin(4);
	t48 = t35 * t37;
	t47 = (t46 + t44 / 0.2e1) * r_i_i_C(1) + (t21 - t22 / 0.2e1) * r_i_i_C(2);
	t26 = pkin(3) + t49;
	t41 = -t42 * r_i_i_C(1) + t10 * r_i_i_C(2) + t53 * t24 + t25 * t26;
	t40 = pkin(2) * cos(t33) + t41;
	t39 = t43 * r_i_i_C(1) + t9 * r_i_i_C(2) - t24 * t26 + t53 * t25;
	t38 = -pkin(2) * sin(t33) + t39;
	t1 = [-sin(qJ(1)) * pkin(1) + t38, t38, t39, (-t24 * t48 - t25 * t36) * pkin(4) + t51, t51; cos(qJ(1)) * pkin(1) + t40, t40, t41, (-t24 * t36 + t25 * t48) * pkin(4) + t52, t52; 0, 0, 0, t34 * t49 + t47, t47;];
	Ja_transl = t1;
end