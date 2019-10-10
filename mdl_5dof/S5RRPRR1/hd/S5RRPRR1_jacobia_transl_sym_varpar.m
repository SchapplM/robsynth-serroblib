% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:59
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRPRR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_jacobia_transl_sym_varpar: pkin has to be [4x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:59:18
	% EndTime: 2019-10-09 20:59:18
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:59:18
	% EndTime: 2019-10-09 20:59:18
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:59:18
	% EndTime: 2019-10-09 20:59:18
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (7->4), mult. (20->10), div. (0->0), fcn. (20->4), ass. (0->7)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t6 = t3 * r_i_i_C(1) - t1 * r_i_i_C(2);
	t5 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t7 = [t4 * r_i_i_C(3) - t6 * t2, t5 * t4, 0, 0, 0; t2 * r_i_i_C(3) + t6 * t4, t5 * t2, 0, 0, 0; 0, t6, 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:59:18
	% EndTime: 2019-10-09 20:59:18
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (15->7), mult. (31->10), div. (0->0), fcn. (33->4), ass. (0->9)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t8 = pkin(1) + r_i_i_C(1);
	t5 = -t1 * r_i_i_C(2) + t8 * t3;
	t7 = r_i_i_C(3) + qJ(3);
	t6 = -r_i_i_C(2) * t3 - t8 * t1;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t5 * t2 + t7 * t4, t6 * t4, t2, 0, 0; t7 * t2 + t5 * t4, t6 * t2, -t4, 0, 0; 0, t5, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:59:18
	% EndTime: 2019-10-09 20:59:18
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (41->10), mult. (41->14), div. (0->0), fcn. (43->6), ass. (0->12)
	t4 = qJ(2) + qJ(4);
	t2 = sin(t4);
	t3 = cos(t4);
	t17 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t10 = pkin(2) + pkin(1);
	t11 = cos(qJ(2)) * t10 + t17;
	t15 = r_i_i_C(3) + pkin(3) + qJ(3);
	t14 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t12 = -t10 * sin(qJ(2)) + t14;
	t9 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [-t11 * t7 + t15 * t9, t12 * t9, t7, t14 * t9, 0; t11 * t9 + t15 * t7, t12 * t7, -t9, t14 * t7, 0; 0, t11, 0, t17, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:59:18
	% EndTime: 2019-10-09 20:59:18
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (88->28), mult. (107->38), div. (0->0), fcn. (117->8), ass. (0->30)
	t16 = qJ(2) + qJ(4);
	t14 = sin(t16);
	t15 = cos(t16);
	t18 = sin(qJ(5));
	t38 = r_i_i_C(2) * t18;
	t44 = pkin(4) + r_i_i_C(3);
	t45 = t14 * t38 + t15 * t44;
	t42 = t44 * t14;
	t24 = pkin(2) + pkin(1);
	t31 = cos(qJ(2)) * t24;
	t41 = t31 + t42;
	t21 = cos(qJ(5));
	t39 = r_i_i_C(1) * t21;
	t20 = sin(qJ(1));
	t35 = t18 * t20;
	t23 = cos(qJ(1));
	t34 = t18 * t23;
	t33 = t20 * t21;
	t32 = t21 * t23;
	t30 = t45 * t20;
	t29 = t14 * t39;
	t27 = t45 * t23;
	t26 = -sin(qJ(2)) * t24 - t29;
	t25 = (-t38 + t39) * t15 + t42;
	t17 = -pkin(3) - qJ(3);
	t4 = t15 * t32 + t35;
	t3 = -t15 * t34 + t33;
	t2 = -t15 * t33 + t34;
	t1 = t15 * t35 + t32;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t23 * t17 - t41 * t20, t26 * t23 + t27, t20, -t23 * t29 + t27, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t20 * t17 + t41 * t23, t26 * t20 + t30, -t23, -t20 * t29 + t30, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, t25 + t31, 0, t25, (-r_i_i_C(1) * t18 - r_i_i_C(2) * t21) * t14;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,5);
end