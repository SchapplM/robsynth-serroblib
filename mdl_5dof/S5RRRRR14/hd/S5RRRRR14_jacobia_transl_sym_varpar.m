% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRR14
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRRR14_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR14_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR14_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
Ja_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:41
	% DurationCPUTime: 0.02s
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
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:41
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (69->17), mult. (74->28), div. (0->0), fcn. (86->8), ass. (0->17)
	t17 = sin(pkin(5));
	t27 = t17 * (pkin(8) + r_i_i_C(3));
	t18 = cos(pkin(5));
	t19 = sin(qJ(3));
	t24 = t18 * t19;
	t20 = cos(qJ(3));
	t23 = t18 * t20;
	t16 = qJ(1) + qJ(2);
	t14 = sin(t16);
	t15 = cos(t16);
	t7 = -t14 * t23 - t15 * t19;
	t8 = -t14 * t24 + t15 * t20;
	t22 = t15 * pkin(2) + t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t14 * t27;
	t5 = t14 * t19 - t15 * t23;
	t6 = -t14 * t20 - t15 * t24;
	t21 = -t14 * pkin(2) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t15 * t27;
	t1 = [-sin(qJ(1)) * pkin(1) + t21, t21, t7 * r_i_i_C(1) - t8 * r_i_i_C(2), 0, 0; cos(qJ(1)) * pkin(1) + t22, t22, -t5 * r_i_i_C(1) + t6 * r_i_i_C(2), 0, 0; 0, 0, (r_i_i_C(1) * t20 - r_i_i_C(2) * t19) * t17, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (221->32), mult. (138->42), div. (0->0), fcn. (134->14), ass. (0->32)
	t31 = sin(pkin(5));
	t32 = cos(pkin(5));
	t33 = sin(qJ(3));
	t48 = -t32 * t33 * pkin(3) + (pkin(8) + pkin(9) + r_i_i_C(3)) * t31;
	t29 = qJ(3) + qJ(4);
	t40 = pkin(5) - t29;
	t39 = sin(t40);
	t24 = pkin(5) + t29;
	t41 = sin(t24) / 0.2e1;
	t13 = t41 - t39 / 0.2e1;
	t30 = qJ(1) + qJ(2);
	t26 = sin(t30);
	t27 = cos(t29);
	t28 = cos(t30);
	t38 = -t13 * t28 - t26 * t27;
	t21 = cos(t24) / 0.2e1;
	t22 = cos(t40);
	t14 = t22 / 0.2e1 + t21;
	t25 = sin(t29);
	t9 = -t14 * t28 + t25 * t26;
	t47 = -t9 * r_i_i_C(1) + t38 * r_i_i_C(2);
	t10 = -t14 * t26 - t25 * t28;
	t37 = t13 * t26 - t27 * t28;
	t46 = t10 * r_i_i_C(1) + t37 * r_i_i_C(2);
	t34 = cos(qJ(3));
	t45 = pkin(3) * t34;
	t43 = t32 * t34;
	t42 = (t41 + t39 / 0.2e1) * r_i_i_C(1) + (t21 - t22 / 0.2e1) * r_i_i_C(2);
	t23 = pkin(2) + t45;
	t36 = -t37 * r_i_i_C(1) + t10 * r_i_i_C(2) + t28 * t23 + t48 * t26;
	t35 = t38 * r_i_i_C(1) + t9 * r_i_i_C(2) - t26 * t23 + t48 * t28;
	t1 = [-sin(qJ(1)) * pkin(1) + t35, t35, (-t26 * t43 - t28 * t33) * pkin(3) + t46, t46, 0; cos(qJ(1)) * pkin(1) + t36, t36, (-t26 * t33 + t28 * t43) * pkin(3) + t47, t47, 0; 0, 0, t31 * t45 + t42, t42, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (386->45), mult. (222->57), div. (0->0), fcn. (216->18), ass. (0->41)
	t41 = cos(qJ(4));
	t30 = t41 * pkin(4) + pkin(3);
	t37 = sin(pkin(5));
	t38 = cos(pkin(5));
	t39 = sin(qJ(4));
	t40 = sin(qJ(3));
	t42 = cos(qJ(3));
	t58 = -(t42 * t39 * pkin(4) + t40 * t30) * t38 + (pkin(8) + pkin(9) + pkin(10) + r_i_i_C(3)) * t37;
	t35 = qJ(3) + qJ(4);
	t34 = qJ(5) + t35;
	t50 = pkin(5) - t34;
	t45 = sin(t50);
	t27 = pkin(5) + t34;
	t51 = sin(t27) / 0.2e1;
	t17 = t51 - t45 / 0.2e1;
	t29 = cos(t34);
	t36 = qJ(1) + qJ(2);
	t32 = sin(t36);
	t33 = cos(t36);
	t49 = -t33 * t17 - t32 * t29;
	t25 = cos(t27) / 0.2e1;
	t26 = cos(t50);
	t18 = t26 / 0.2e1 + t25;
	t28 = sin(t34);
	t9 = -t33 * t18 + t32 * t28;
	t57 = -t9 * r_i_i_C(1) + t49 * r_i_i_C(2);
	t10 = -t32 * t18 - t33 * t28;
	t48 = t32 * t17 - t33 * t29;
	t56 = t10 * r_i_i_C(1) + t48 * r_i_i_C(2);
	t55 = pkin(4) * sin(t35);
	t53 = t39 * t40;
	t52 = (t51 + t45 / 0.2e1) * r_i_i_C(1) + (t25 - t26 / 0.2e1) * r_i_i_C(2);
	t47 = -pkin(4) * t53 + t30 * t42;
	t19 = pkin(2) + pkin(4) * cos(t35) + t42 * pkin(3);
	t46 = -t48 * r_i_i_C(1) + t10 * r_i_i_C(2) + t33 * t19 + t58 * t32;
	t44 = pkin(4) * (t41 * t42 - t53);
	t43 = t49 * r_i_i_C(1) + t9 * r_i_i_C(2) - t32 * t19 + t58 * t33;
	t20 = -t40 * pkin(3) - t55;
	t13 = t38 * t44;
	t12 = t47 * t38;
	t1 = [-sin(qJ(1)) * pkin(1) + t43, t43, -t32 * t12 + t33 * t20 + t56, -t32 * t13 - t33 * t55 + t56, t56; cos(qJ(1)) * pkin(1) + t46, t46, t33 * t12 + t32 * t20 + t57, t33 * t13 - t32 * t55 + t57, t57; 0, 0, t47 * t37 + t52, t37 * t44 + t52, t52;];
	Ja_transl = t1;
end