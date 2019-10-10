% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPP3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
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
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (11->6), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->9)
	t8 = pkin(7) + r_i_i_C(3);
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t2 * t5 + t4 * t8, t6 * t4, 0, 0, 0, 0; t2 * t8 + t4 * t5, t6 * t2, 0, 0, 0, 0; 0, t7, 0, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (37->9), mult. (41->14), div. (0->0), fcn. (41->6), ass. (0->12)
	t6 = qJ(2) + qJ(3);
	t3 = sin(t6);
	t4 = cos(t6);
	t14 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t17 = t14 + cos(qJ(2)) * pkin(2);
	t15 = r_i_i_C(3) + pkin(8) + pkin(7);
	t13 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t12 = pkin(1) + t17;
	t11 = -sin(qJ(2)) * pkin(2) + t13;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t12 * t8 + t15 * t9, t11 * t9, t13 * t9, 0, 0, 0; t12 * t9 + t15 * t8, t11 * t8, t13 * t8, 0, 0, 0; 0, t17, t14, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (100->28), mult. (121->39), div. (0->0), fcn. (129->8), ass. (0->29)
	t19 = qJ(2) + qJ(3);
	t16 = sin(t19);
	t17 = cos(t19);
	t20 = sin(qJ(4));
	t39 = r_i_i_C(2) * t20;
	t45 = pkin(9) + r_i_i_C(3);
	t46 = t16 * t39 + t17 * t45;
	t43 = t17 * pkin(3) + t45 * t16;
	t18 = cos(qJ(2)) * pkin(2);
	t42 = pkin(1) + t18 + t43;
	t23 = cos(qJ(4));
	t40 = r_i_i_C(1) * t23;
	t22 = sin(qJ(1));
	t36 = t20 * t22;
	t24 = cos(qJ(1));
	t35 = t20 * t24;
	t34 = t22 * t23;
	t33 = t23 * t24;
	t32 = t46 * t22;
	t30 = t46 * t24;
	t28 = (-pkin(3) - t40) * t16;
	t27 = (-t39 + t40) * t17 + t43;
	t26 = -sin(qJ(2)) * pkin(2) + t28;
	t25 = -pkin(8) - pkin(7);
	t4 = t17 * t33 + t36;
	t3 = -t17 * t35 + t34;
	t2 = -t17 * t34 + t35;
	t1 = t17 * t36 + t33;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t42 * t22 - t24 * t25, t26 * t24 + t30, t24 * t28 + t30, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t22 * t25 + t42 * t24, t26 * t22 + t32, t22 * t28 + t32, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, t18 + t27, t27, (-r_i_i_C(1) * t20 - r_i_i_C(2) * t23) * t16, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (157->31), mult. (196->42), div. (0->0), fcn. (217->8), ass. (0->30)
	t23 = qJ(2) + qJ(3);
	t20 = sin(t23);
	t21 = cos(t23);
	t27 = cos(qJ(4));
	t52 = pkin(9) + r_i_i_C(1);
	t53 = r_i_i_C(2) * t20 * t27 + t21 * t52;
	t47 = pkin(4) - r_i_i_C(2);
	t24 = sin(qJ(4));
	t37 = r_i_i_C(3) + qJ(5);
	t51 = t37 * t24;
	t49 = t21 * pkin(3) + t52 * t20;
	t22 = cos(qJ(2)) * pkin(2);
	t48 = pkin(1) + t22 + t49;
	t26 = sin(qJ(1));
	t41 = t24 * t26;
	t40 = t26 * t27;
	t28 = cos(qJ(1));
	t39 = t27 * t28;
	t38 = t28 * t24;
	t36 = t53 * t26;
	t34 = t53 * t28;
	t32 = (-pkin(4) * t27 - pkin(3) - t51) * t20;
	t31 = t49 + (t27 * t47 + t51) * t21;
	t30 = -sin(qJ(2)) * pkin(2) + t32;
	t29 = -pkin(8) - pkin(7);
	t4 = t21 * t39 + t41;
	t3 = t21 * t38 - t40;
	t2 = t21 * t40 - t38;
	t1 = t21 * t41 + t39;
	t5 = [-t37 * t1 - t47 * t2 - t48 * t26 - t28 * t29, t30 * t28 + t34, t28 * t32 + t34, -t47 * t3 + t37 * t4, t3, 0; -t26 * t29 + t48 * t28 + t37 * t3 + t47 * t4, t30 * t26 + t36, t26 * t32 + t36, -t47 * t1 + t37 * t2, t1, 0; 0, t22 + t31, t31, (-t47 * t24 + t37 * t27) * t20, t20 * t24, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (204->29), mult. (251->40), div. (0->0), fcn. (281->8), ass. (0->30)
	t27 = sin(qJ(4));
	t30 = cos(qJ(4));
	t38 = pkin(4) + r_i_i_C(3) + qJ(6);
	t41 = r_i_i_C(2) + qJ(5);
	t56 = t41 * t27 + t38 * t30;
	t53 = pkin(5) + pkin(9) + r_i_i_C(1);
	t26 = qJ(2) + qJ(3);
	t24 = cos(t26);
	t52 = t24 * t53;
	t23 = sin(t26);
	t51 = t24 * pkin(3) + t53 * t23;
	t25 = cos(qJ(2)) * pkin(2);
	t50 = pkin(1) + t25 + t51;
	t29 = sin(qJ(1));
	t45 = t29 * t27;
	t44 = t29 * t30;
	t31 = cos(qJ(1));
	t43 = t31 * t27;
	t42 = t31 * t30;
	t39 = t29 * t52;
	t37 = t31 * t52;
	t35 = t56 * t24 + t51;
	t34 = (-pkin(3) - t56) * t23;
	t33 = -sin(qJ(2)) * pkin(2) + t34;
	t32 = -pkin(8) - pkin(7);
	t4 = t24 * t42 + t45;
	t3 = t24 * t43 - t44;
	t2 = t24 * t44 - t43;
	t1 = t24 * t45 + t42;
	t5 = [-t41 * t1 - t38 * t2 - t50 * t29 - t31 * t32, t31 * t33 + t37, t31 * t34 + t37, -t38 * t3 + t41 * t4, t3, t4; -t29 * t32 + t41 * t3 + t50 * t31 + t38 * t4, t29 * t33 + t39, t29 * t34 + t39, -t38 * t1 + t41 * t2, t1, t2; 0, t25 + t35, t35, (-t38 * t27 + t41 * t30) * t23, t23 * t27, t23 * t30;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end