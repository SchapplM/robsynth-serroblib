% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:02
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRRR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_jacobia_transl_sym_varpar: pkin has to be [2x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:47
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:47
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
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:47
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:47
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (29->7), mult. (32->12), div. (0->0), fcn. (32->6), ass. (0->10)
	t8 = sin(qJ(3));
	t9 = cos(qJ(3));
	t15 = t9 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t12 = -r_i_i_C(1) * t8 - r_i_i_C(2) * t9;
	t7 = qJ(1) + qJ(2);
	t5 = sin(t7);
	t6 = cos(t7);
	t11 = t6 * r_i_i_C(3) - t15 * t5;
	t10 = t5 * r_i_i_C(3) + t15 * t6;
	t1 = [-sin(qJ(1)) * pkin(1) + t11, t11, t12 * t6, 0, 0; cos(qJ(1)) * pkin(1) + t10, t10, t12 * t5, 0, 0; 0, 0, t15, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:47
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (69->10), mult. (55->16), div. (0->0), fcn. (55->8), ass. (0->13)
	t11 = qJ(3) + qJ(4);
	t7 = sin(t11);
	t9 = cos(t11);
	t19 = t9 * r_i_i_C(1) - t7 * r_i_i_C(2);
	t23 = t19 + cos(qJ(3)) * pkin(2);
	t18 = -r_i_i_C(1) * t7 - r_i_i_C(2) * t9;
	t12 = qJ(1) + qJ(2);
	t10 = cos(t12);
	t8 = sin(t12);
	t17 = t8 * r_i_i_C(3) + t10 * t23;
	t16 = -sin(qJ(3)) * pkin(2) + t18;
	t15 = t10 * r_i_i_C(3) - t23 * t8;
	t1 = [-sin(qJ(1)) * pkin(1) + t15, t15, t16 * t10, t18 * t10, 0; cos(qJ(1)) * pkin(1) + t17, t17, t16 * t8, t18 * t8, 0; 0, 0, t23, t19, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:47
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (126->25), mult. (117->40), div. (0->0), fcn. (129->10), ass. (0->28)
	t21 = qJ(3) + qJ(4);
	t17 = sin(t21);
	t16 = t17 * r_i_i_C(3);
	t40 = cos(qJ(3)) * pkin(2);
	t42 = t16 + t40;
	t19 = cos(t21);
	t23 = sin(qJ(5));
	t41 = r_i_i_C(2) * t17 * t23 + r_i_i_C(3) * t19;
	t22 = qJ(1) + qJ(2);
	t18 = sin(t22);
	t38 = t41 * t18;
	t37 = t19 * t23;
	t25 = cos(qJ(5));
	t36 = t19 * t25;
	t20 = cos(t22);
	t35 = t20 * t23;
	t34 = t20 * t25;
	t33 = t41 * t20;
	t32 = r_i_i_C(1) * t17 * t25;
	t7 = t18 * t25 - t19 * t35;
	t8 = t18 * t23 + t19 * t34;
	t30 = t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t42 * t20;
	t29 = r_i_i_C(1) * t36 - r_i_i_C(2) * t37 + t16;
	t28 = -sin(qJ(3)) * pkin(2) - t32;
	t5 = t18 * t37 + t34;
	t6 = -t18 * t36 + t35;
	t27 = t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t42 * t18;
	t1 = [-sin(qJ(1)) * pkin(1) + t27, t27, t28 * t20 + t33, -t20 * t32 + t33, r_i_i_C(1) * t7 - r_i_i_C(2) * t8; cos(qJ(1)) * pkin(1) + t30, t30, t28 * t18 + t38, -t18 * t32 + t38, -r_i_i_C(1) * t5 + r_i_i_C(2) * t6; 0, 0, t29 + t40, t29, (-r_i_i_C(1) * t23 - r_i_i_C(2) * t25) * t17;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end