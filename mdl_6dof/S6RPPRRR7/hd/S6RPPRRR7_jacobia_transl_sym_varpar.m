% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRR7
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:11
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRR7_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR7_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR7_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (7->5), mult. (8->4), div. (0->0), fcn. (10->2), ass. (0->5)
	t4 = pkin(1) - r_i_i_C(2);
	t3 = r_i_i_C(3) + qJ(2);
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t5 = [-t1 * t4 + t2 * t3, t1, 0, 0, 0, 0; t1 * t3 + t2 * t4, -t2, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (11->6), mult. (16->6), div. (0->0), fcn. (20->4), ass. (0->5)
	t6 = pkin(1) + r_i_i_C(3) + qJ(3);
	t5 = r_i_i_C(1) * sin(pkin(10)) + r_i_i_C(2) * cos(pkin(10)) + qJ(2);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t6 * t3 + t5 * t4, t3, t4, 0, 0, 0; t5 * t3 + t6 * t4, -t4, t3, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (28->10), mult. (30->11), div. (0->0), fcn. (34->5), ass. (0->10)
	t11 = pkin(1) + r_i_i_C(3) + pkin(7) + qJ(3);
	t3 = pkin(10) + qJ(4);
	t1 = sin(t3);
	t2 = cos(t3);
	t10 = r_i_i_C(1) * t2 - r_i_i_C(2) * t1;
	t9 = -t1 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t8 = sin(pkin(10)) * pkin(3) + qJ(2) - t9;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t4 = [-t11 * t6 + t8 * t7, t6, t7, t10 * t6, 0, 0; t11 * t7 + t8 * t6, -t7, t6, -t10 * t7, 0, 0; 0, 0, 0, t9, 0, 0;];
	Ja_transl = t4;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (65->16), mult. (47->17), div. (0->0), fcn. (51->7), ass. (0->16)
	t10 = pkin(10) + qJ(4);
	t8 = qJ(5) + t10;
	t4 = sin(t8);
	t5 = cos(t8);
	t14 = -t4 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t20 = t14 - pkin(4) * sin(t10);
	t18 = pkin(4) * cos(t10);
	t17 = r_i_i_C(1) * t5;
	t16 = r_i_i_C(2) * t4;
	t15 = pkin(1) + r_i_i_C(3) + pkin(8) + pkin(7) + qJ(3);
	t13 = qJ(2) + sin(pkin(10)) * pkin(3) - t20;
	t12 = cos(qJ(1));
	t11 = sin(qJ(1));
	t3 = t12 * t16;
	t2 = t11 * t17;
	t1 = [-t15 * t11 + t13 * t12, t11, t12, t2 + (-t16 + t18) * t11, -t11 * t16 + t2, 0; t13 * t11 + t15 * t12, -t12, t11, t3 + (-t17 - t18) * t12, -t12 * t17 + t3, 0; 0, 0, 0, t20, t14, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (158->35), mult. (127->43), div. (0->0), fcn. (139->9), ass. (0->31)
	t20 = pkin(10) + qJ(4);
	t18 = qJ(5) + t20;
	t14 = sin(t18);
	t41 = pkin(9) + r_i_i_C(3);
	t43 = t41 * t14;
	t15 = cos(t18);
	t42 = t41 * t15;
	t39 = pkin(4) * sin(t20);
	t38 = pkin(4) * cos(t20);
	t21 = sin(qJ(6));
	t37 = r_i_i_C(2) * t21;
	t36 = pkin(1) + pkin(8) + pkin(7) + qJ(3);
	t24 = cos(qJ(1));
	t34 = t21 * t24;
	t22 = sin(qJ(1));
	t33 = t22 * t21;
	t23 = cos(qJ(6));
	t32 = t22 * t23;
	t31 = t23 * t24;
	t30 = t15 * t37;
	t29 = -r_i_i_C(1) * t23 - pkin(5);
	t28 = t22 * t43 + (pkin(5) * t22 + r_i_i_C(1) * t32) * t15;
	t27 = t42 + (t29 + t37) * t14;
	t26 = pkin(5) * t14 - t42 + qJ(2) + t39 + sin(pkin(10)) * pkin(3);
	t25 = t29 * t15 - t43;
	t6 = t24 * t30;
	t4 = t14 * t31 - t33;
	t3 = t14 * t34 + t32;
	t2 = t14 * t32 + t34;
	t1 = -t14 * t33 + t31;
	t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t36 * t22 + t26 * t24, t22, t24, (-t30 + t38) * t22 + t28, -t22 * t30 + t28, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t26 * t22 + t36 * t24, -t24, t22, t6 + (t25 - t38) * t24, t25 * t24 + t6, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 0, 0, t27 - t39, t27, (-r_i_i_C(1) * t21 - r_i_i_C(2) * t23) * t15;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end