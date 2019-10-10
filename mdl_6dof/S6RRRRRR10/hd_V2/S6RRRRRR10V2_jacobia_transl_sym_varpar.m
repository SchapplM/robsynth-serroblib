% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:38
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR10V2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10V2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR10V2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobia_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:01
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:01
	% EndTime: 2019-10-10 13:38:02
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
	% StartTime: 2019-10-10 13:38:01
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (9->5), mult. (22->10), div. (0->0), fcn. (22->4), ass. (0->8)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t8 = [t4 * r_i_i_C(3) - t2 * t5, t6 * t4, 0, 0, 0, 0; t2 * r_i_i_C(3) + t4 * t5, t6 * t2, 0, 0, 0, 0; 0, t7, 0, 0, 0, 0;];
	Ja_transl = t8;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:01
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (33->8), mult. (39->14), div. (0->0), fcn. (39->6), ass. (0->11)
	t6 = qJ(2) + qJ(3);
	t3 = sin(t6);
	t4 = cos(t6);
	t13 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t15 = t13 + cos(qJ(2)) * pkin(2);
	t12 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t11 = pkin(1) + t15;
	t10 = -sin(qJ(2)) * pkin(2) + t12;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [t9 * r_i_i_C(3) - t11 * t8, t10 * t9, t12 * t9, 0, 0, 0; t8 * r_i_i_C(3) + t11 * t9, t10 * t8, t12 * t8, 0, 0, 0; 0, t15, t13, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:01
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (96->25), mult. (119->37), div. (0->0), fcn. (127->8), ass. (0->28)
	t19 = qJ(2) + qJ(3);
	t16 = sin(t19);
	t17 = cos(t19);
	t20 = sin(qJ(4));
	t38 = r_i_i_C(2) * t20;
	t44 = pkin(5) + r_i_i_C(3);
	t45 = t16 * t38 + t17 * t44;
	t42 = t17 * pkin(3) + t44 * t16;
	t18 = cos(qJ(2)) * pkin(2);
	t41 = pkin(1) + t18 + t42;
	t23 = cos(qJ(4));
	t39 = r_i_i_C(1) * t23;
	t24 = cos(qJ(1));
	t35 = t20 * t24;
	t22 = sin(qJ(1));
	t34 = t22 * t20;
	t33 = t22 * t23;
	t32 = t23 * t24;
	t31 = t45 * t22;
	t29 = t45 * t24;
	t27 = (-pkin(3) - t39) * t16;
	t26 = (-t38 + t39) * t17 + t42;
	t25 = -sin(qJ(2)) * pkin(2) + t27;
	t4 = t17 * t32 + t34;
	t3 = -t17 * t35 + t33;
	t2 = -t17 * t33 + t35;
	t1 = t17 * t34 + t32;
	t5 = [r_i_i_C(1) * t2 + r_i_i_C(2) * t1 - t41 * t22, t25 * t24 + t29, t24 * t27 + t29, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t41 * t24, t25 * t22 + t31, t22 * t27 + t31, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, t18 + t26, t26, (-r_i_i_C(1) * t20 - r_i_i_C(2) * t23) * t16, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:01
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (163->43), mult. (225->70), div. (0->0), fcn. (263->10), ass. (0->39)
	t24 = qJ(2) + qJ(3);
	t22 = cos(t24);
	t21 = sin(t24);
	t29 = cos(qJ(5));
	t25 = sin(qJ(5));
	t30 = cos(qJ(4));
	t48 = t25 * t30;
	t35 = t21 * t48 + t22 * t29;
	t45 = t29 * t30;
	t36 = -t21 * t45 + t22 * t25;
	t54 = pkin(5) * t22 + r_i_i_C(1) * t36 + r_i_i_C(2) * t35;
	t26 = sin(qJ(4));
	t52 = r_i_i_C(3) * t26;
	t16 = t21 * pkin(5);
	t17 = t22 * pkin(3);
	t51 = t21 * t25;
	t50 = t21 * t29;
	t31 = cos(qJ(1));
	t49 = t21 * t31;
	t28 = sin(qJ(1));
	t47 = t28 * t26;
	t46 = t28 * t30;
	t44 = t31 * t26;
	t43 = t31 * t30;
	t42 = t54 * t28;
	t41 = t54 * t31;
	t23 = cos(qJ(2)) * pkin(2);
	t40 = -t23 - pkin(1) - t17;
	t39 = (t22 * t45 + t51) * r_i_i_C(1) + (-t22 * t48 + t50) * r_i_i_C(2) + t22 * t52 + t17 + t16;
	t38 = (-pkin(3) - t52) * t21;
	t37 = t29 * r_i_i_C(1) - t25 * r_i_i_C(2);
	t32 = -sin(qJ(2)) * pkin(2) + t38;
	t12 = t22 * t43 + t47;
	t11 = t22 * t44 - t46;
	t10 = t22 * t46 - t44;
	t9 = -t22 * t47 - t43;
	t2 = t12 * t29 + t25 * t49;
	t1 = -t12 * t25 + t29 * t49;
	t3 = [t9 * r_i_i_C(3) - t37 * t10 + ((-t25 * r_i_i_C(1) - t29 * r_i_i_C(2) - pkin(5)) * t21 + t40) * t28, t32 * t31 + t41, t31 * t38 + t41, t12 * r_i_i_C(3) - t37 * t11, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t11 * r_i_i_C(3) + (-t40 + t16) * t31, t32 * t28 + t42, t28 * t38 + t42, t10 * r_i_i_C(3) + t37 * t9, (-t10 * t25 + t28 * t50) * r_i_i_C(1) + (-t10 * t29 - t28 * t51) * r_i_i_C(2), 0; 0, t23 + t39, t39, (r_i_i_C(3) * t30 - t37 * t26) * t21, -t35 * r_i_i_C(1) + t36 * r_i_i_C(2), 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:01
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (342->62), mult. (503->112), div. (0->0), fcn. (621->12), ass. (0->56)
	t42 = qJ(2) + qJ(3);
	t40 = cos(t42);
	t39 = sin(t42);
	t49 = cos(qJ(5));
	t44 = sin(qJ(5));
	t50 = cos(qJ(4));
	t73 = t44 * t50;
	t54 = t39 * t73 + t40 * t49;
	t81 = r_i_i_C(3) + pkin(6);
	t86 = pkin(5) * t40 - t54 * t81;
	t43 = sin(qJ(6));
	t48 = cos(qJ(6));
	t60 = t48 * r_i_i_C(1) - t43 * r_i_i_C(2);
	t45 = sin(qJ(4));
	t51 = cos(qJ(1));
	t67 = t51 * t45;
	t47 = sin(qJ(1));
	t69 = t47 * t50;
	t27 = t40 * t69 - t67;
	t77 = t39 * t47;
	t10 = t27 * t49 + t44 * t77;
	t66 = t51 * t50;
	t70 = t47 * t45;
	t26 = t40 * t70 + t66;
	t84 = t10 * t43 - t26 * t48;
	t83 = -t10 * t48 - t26 * t43;
	t82 = t40 * pkin(3) + t39 * pkin(5);
	t76 = t39 * t49;
	t75 = t39 * t51;
	t74 = t43 * t45;
	t72 = t45 * t48;
	t71 = t45 * t49;
	t68 = t49 * t50;
	t65 = t39 * t74;
	t64 = t39 * t72;
	t63 = t39 * t67;
	t62 = t81 * t44;
	t61 = -sin(qJ(2)) * pkin(2) - pkin(3) * t39;
	t59 = t43 * r_i_i_C(1) + t48 * r_i_i_C(2);
	t41 = cos(qJ(2)) * pkin(2);
	t58 = t41 + pkin(1) + t82;
	t23 = t39 * t68 - t40 * t44;
	t57 = (-t65 * r_i_i_C(1) - t64 * r_i_i_C(2) - t60 * t23 + t86) * t47;
	t20 = t23 * t51;
	t56 = (t20 * t43 - t48 * t63) * r_i_i_C(2) + (-t20 * t48 - t43 * t63) * r_i_i_C(1) + t86 * t51;
	t55 = -t27 * t44 + t47 * t76;
	t25 = t39 * t44 + t40 * t68;
	t53 = (-t25 * t43 + t40 * t72) * r_i_i_C(2) + (t25 * t48 + t40 * t74) * r_i_i_C(1) + t82 + t81 * (t40 * t73 - t76);
	t52 = -t60 * t49 - t62;
	t29 = t40 * t66 + t70;
	t28 = t40 * t67 - t69;
	t14 = t29 * t49 + t44 * t75;
	t13 = t29 * t44 - t49 * t75;
	t2 = t14 * t48 + t28 * t43;
	t1 = -t14 * t43 + t28 * t48;
	t3 = [t83 * r_i_i_C(1) + t84 * r_i_i_C(2) - t58 * t47 + t81 * t55, t61 * t51 + t56, -pkin(3) * t75 + t56, t52 * t28 + t59 * t29, -t60 * t13 + t81 * t14, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t81 * t13 + t58 * t51, t61 * t47 + t57, -pkin(3) * t77 + t57, t52 * t26 + t59 * t27, t81 * t10 + t60 * t55, -t84 * r_i_i_C(1) + t83 * r_i_i_C(2); 0, t41 + t53, t53, ((t43 * t50 - t48 * t71) * r_i_i_C(1) + (t43 * t71 + t48 * t50) * r_i_i_C(2) - t45 * t62) * t39, t81 * t23 - t60 * t54, (-t23 * t43 + t64) * r_i_i_C(1) + (-t23 * t48 - t65) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end