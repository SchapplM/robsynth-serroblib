% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRR11
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
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRRR11_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR11_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR11_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
Ja_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	% Symbolic code from jacobia_transl_1_floatb_twist_matlab.m not found
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (6->3), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->4)
	t3 = pkin(10) + qJ(2);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [0, -r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; 0, r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (41->15), mult. (46->24), div. (0->0), fcn. (54->6), ass. (0->17)
	t10 = sin(qJ(3));
	t7 = pkin(10) + qJ(2);
	t5 = sin(t7);
	t16 = t5 * t10;
	t11 = cos(qJ(3));
	t15 = t5 * t11;
	t6 = cos(t7);
	t14 = t6 * t10;
	t13 = t6 * t11;
	t8 = sin(pkin(5));
	t12 = (pkin(7) + r_i_i_C(3)) * t8;
	t9 = cos(pkin(5));
	t4 = -t9 * t16 + t13;
	t3 = -t9 * t15 - t14;
	t2 = -t9 * t14 - t15;
	t1 = -t9 * t13 + t16;
	t17 = [0, -t5 * pkin(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t6 * t12, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0; 0, t6 * pkin(2) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t5 * t12, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0; 1, 0, (r_i_i_C(1) * t11 - r_i_i_C(2) * t10) * t8, 0, 0;];
	Ja_transl = t17;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (163->30), mult. (100->40), div. (0->0), fcn. (96->12), ass. (0->28)
	t22 = pkin(10) + qJ(2);
	t17 = sin(t22);
	t18 = cos(t22);
	t23 = qJ(3) + qJ(4);
	t21 = cos(t23);
	t19 = pkin(5) + t23;
	t14 = cos(t19) / 0.2e1;
	t30 = pkin(5) - t23;
	t15 = cos(t30);
	t10 = t15 / 0.2e1 + t14;
	t20 = sin(t23);
	t5 = -t10 * t18 + t17 * t20;
	t29 = sin(t30);
	t32 = sin(t19) / 0.2e1;
	t9 = t32 - t29 / 0.2e1;
	t37 = -t5 * r_i_i_C(1) + (-t17 * t21 - t18 * t9) * r_i_i_C(2);
	t6 = -t10 * t17 - t18 * t20;
	t36 = t6 * r_i_i_C(1) + (t17 * t9 - t18 * t21) * r_i_i_C(2);
	t35 = (t32 + t29 / 0.2e1) * r_i_i_C(1) + (t14 - t15 / 0.2e1) * r_i_i_C(2);
	t27 = cos(qJ(3));
	t34 = pkin(3) * t27;
	t25 = cos(pkin(5));
	t33 = t25 * t27;
	t31 = r_i_i_C(1) * t21 + pkin(2) + t34;
	t24 = sin(pkin(5));
	t26 = sin(qJ(3));
	t28 = -t25 * t26 * pkin(3) - r_i_i_C(1) * t9 + (r_i_i_C(3) + pkin(7) + pkin(8)) * t24;
	t1 = [0, t5 * r_i_i_C(2) - t31 * t17 + t28 * t18, (-t17 * t33 - t18 * t26) * pkin(3) + t36, t36, 0; 0, t6 * r_i_i_C(2) + t28 * t17 + t31 * t18, (-t17 * t26 + t18 * t33) * pkin(3) + t37, t37, 0; 1, 0, t24 * t34 + t35, t35, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (306->43), mult. (174->55), div. (0->0), fcn. (168->16), ass. (0->37)
	t29 = qJ(3) + qJ(4);
	t27 = qJ(5) + t29;
	t40 = pkin(5) - t27;
	t37 = sin(t40);
	t20 = pkin(5) + t27;
	t42 = sin(t20) / 0.2e1;
	t12 = t42 - t37 / 0.2e1;
	t22 = cos(t27);
	t28 = pkin(10) + qJ(2);
	t24 = sin(t28);
	t25 = cos(t28);
	t18 = cos(t20) / 0.2e1;
	t19 = cos(t40);
	t13 = t19 / 0.2e1 + t18;
	t21 = sin(t27);
	t5 = -t25 * t13 + t24 * t21;
	t47 = -t5 * r_i_i_C(1) + (-t25 * t12 - t24 * t22) * r_i_i_C(2);
	t6 = -t24 * t13 - t25 * t21;
	t46 = t6 * r_i_i_C(1) + (t24 * t12 - t25 * t22) * r_i_i_C(2);
	t45 = pkin(4) * sin(t29);
	t32 = sin(qJ(4));
	t33 = sin(qJ(3));
	t44 = t32 * t33;
	t43 = (t42 + t37 / 0.2e1) * r_i_i_C(1) + (t18 - t19 / 0.2e1) * r_i_i_C(2);
	t35 = cos(qJ(3));
	t41 = t22 * r_i_i_C(1) + pkin(2) + pkin(4) * cos(t29) + t35 * pkin(3);
	t34 = cos(qJ(4));
	t23 = t34 * pkin(4) + pkin(3);
	t30 = sin(pkin(5));
	t31 = cos(pkin(5));
	t39 = -t12 * r_i_i_C(1) - (t35 * t32 * pkin(4) + t33 * t23) * t31 + (r_i_i_C(3) + pkin(8) + pkin(9) + pkin(7)) * t30;
	t38 = -pkin(4) * t44 + t23 * t35;
	t36 = pkin(4) * (t34 * t35 - t44);
	t15 = -t33 * pkin(3) - t45;
	t9 = t31 * t36;
	t8 = t38 * t31;
	t1 = [0, t5 * r_i_i_C(2) - t41 * t24 + t39 * t25, t25 * t15 - t24 * t8 + t46, -t24 * t9 - t25 * t45 + t46, t46; 0, t6 * r_i_i_C(2) + t39 * t24 + t41 * t25, t24 * t15 + t25 * t8 + t47, -t24 * t45 + t25 * t9 + t47, t47; 1, 0, t38 * t30 + t43, t30 * t36 + t43, t43;];
	Ja_transl = t1;
end