% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRR3
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRRR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobia_transl_sym_varpar: pkin has to be [5x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:04:40
	% EndTime: 2019-10-09 21:04:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:04:40
	% EndTime: 2019-10-09 21:04:40
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
	% StartTime: 2019-10-09 21:04:40
	% EndTime: 2019-10-09 21:04:40
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
	% StartTime: 2019-10-09 21:04:40
	% EndTime: 2019-10-09 21:04:40
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (31->7), mult. (39->14), div. (0->0), fcn. (39->6), ass. (0->10)
	t4 = qJ(2) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t14 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t9 = cos(qJ(2)) * pkin(1) + t14;
	t11 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t10 = -sin(qJ(2)) * pkin(1) + t11;
	t8 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [t8 * r_i_i_C(3) - t9 * t6, t10 * t8, t11 * t8, 0, 0; t6 * r_i_i_C(3) + t9 * t8, t10 * t6, t11 * t6, 0, 0; 0, t9, t14, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:04:40
	% EndTime: 2019-10-09 21:04:40
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (94->25), mult. (119->37), div. (0->0), fcn. (127->8), ass. (0->28)
	t17 = qJ(2) + qJ(3);
	t15 = sin(t17);
	t16 = cos(t17);
	t18 = sin(qJ(4));
	t37 = r_i_i_C(2) * t18;
	t44 = pkin(5) + r_i_i_C(3);
	t45 = t15 * t37 + t16 * t44;
	t42 = t16 * pkin(2) + t44 * t15;
	t39 = cos(qJ(2)) * pkin(1);
	t41 = t39 + t42;
	t21 = cos(qJ(4));
	t38 = r_i_i_C(1) * t21;
	t23 = cos(qJ(1));
	t34 = t18 * t23;
	t20 = sin(qJ(1));
	t33 = t20 * t18;
	t32 = t20 * t21;
	t31 = t21 * t23;
	t30 = t45 * t20;
	t28 = t45 * t23;
	t27 = (-pkin(2) - t38) * t15;
	t25 = (-t37 + t38) * t16 + t42;
	t24 = -sin(qJ(2)) * pkin(1) + t27;
	t4 = t16 * t31 + t33;
	t3 = -t16 * t34 + t32;
	t2 = -t16 * t32 + t34;
	t1 = t16 * t33 + t31;
	t5 = [r_i_i_C(1) * t2 + r_i_i_C(2) * t1 - t41 * t20, t23 * t24 + t28, t23 * t27 + t28, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t41 * t23, t20 * t24 + t30, t20 * t27 + t30, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, t25 + t39, t25, (-r_i_i_C(1) * t18 - r_i_i_C(2) * t21) * t15, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:04:40
	% EndTime: 2019-10-09 21:04:40
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (165->35), mult. (163->52), div. (0->0), fcn. (175->10), ass. (0->36)
	t54 = pkin(5) + r_i_i_C(3);
	t29 = cos(qJ(4));
	t19 = pkin(3) * t29 + pkin(2);
	t25 = qJ(2) + qJ(3);
	t21 = sin(t25);
	t23 = cos(t25);
	t53 = t23 * t19 + t54 * t21;
	t48 = cos(qJ(2)) * pkin(1);
	t52 = t48 + t53;
	t24 = qJ(4) + qJ(5);
	t22 = cos(t24);
	t31 = cos(qJ(1));
	t20 = sin(t24);
	t28 = sin(qJ(1));
	t41 = t28 * t20;
	t5 = t22 * t31 + t23 * t41;
	t40 = t28 * t22;
	t6 = t20 * t31 - t23 * t40;
	t50 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t44 = t23 * t31;
	t7 = -t20 * t44 + t40;
	t8 = t22 * t44 + t41;
	t49 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t47 = r_i_i_C(1) * t22;
	t46 = r_i_i_C(2) * t20;
	t26 = sin(qJ(4));
	t43 = t26 * t28;
	t42 = t26 * t31;
	t39 = t21 * t46;
	t38 = (t23 * t54 + t39) * t28;
	t37 = t31 * t39 + t54 * t44;
	t36 = -r_i_i_C(1) * t20 - r_i_i_C(2) * t22;
	t35 = (-t19 - t47) * t21;
	t33 = (-t46 + t47) * t23 + t53;
	t32 = -sin(qJ(2)) * pkin(1) + t35;
	t1 = [pkin(3) * t42 + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t52 * t28, t32 * t31 + t37, t31 * t35 + t37, (-t23 * t42 + t28 * t29) * pkin(3) + t49, t49; pkin(3) * t43 + t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t52 * t31, t32 * t28 + t38, t28 * t35 + t38, (-t23 * t43 - t29 * t31) * pkin(3) + t50, t50; 0, t33 + t48, t33, (-pkin(3) * t26 + t36) * t21, t36 * t21;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end