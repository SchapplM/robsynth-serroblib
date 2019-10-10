% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6PRPPRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:28
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPPRR3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (69->16), mult. (180->37), div. (39->10), fcn. (276->9), ass. (0->23)
	t36 = sin(pkin(6));
	t40 = cos(qJ(2));
	t44 = t36 * t40;
	t38 = cos(pkin(6));
	t39 = sin(qJ(2));
	t43 = t38 * t39;
	t42 = t38 * t40;
	t41 = t36 ^ 2;
	t37 = cos(pkin(10));
	t35 = sin(pkin(10));
	t34 = 0.1e1 / t40 ^ 2;
	t31 = -t35 * t43 + t37 * t40;
	t30 = t35 * t42 + t37 * t39;
	t29 = t35 * t40 + t37 * t43;
	t27 = t35 * t39 - t37 * t42;
	t26 = 0.1e1 / t31 ^ 2;
	t24 = atan2(-t27, -t44);
	t23 = cos(t24);
	t22 = sin(t24);
	t21 = -t22 * t27 - t23 * t44;
	t20 = 0.1e1 / t21 ^ 2;
	t18 = (t29 / t40 + t39 * t27 * t34) / t36 / (0.1e1 + t27 ^ 2 / t41 * t34);
	t1 = [0, t18, 0, 0, 0, 0; 0, (t31 / t21 - (t23 * t36 * t39 - t22 * t29 + (t22 * t44 - t23 * t27) * t18) * t30 * t20) / (t30 ^ 2 * t20 + 0.1e1), 0, 0, 0, 0; 0, -t30 * t35 * t36 * t26 / (t35 ^ 2 * t41 * t26 + 0.1e1), 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->0), mult. (78->0), div. (6->0), fcn. (108->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, -1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (222->23), mult. (626->54), div. (35->9), fcn. (880->13), ass. (0->36)
	t60 = sin(pkin(10));
	t61 = sin(pkin(6));
	t73 = t60 * t61;
	t64 = cos(pkin(6));
	t66 = sin(qJ(2));
	t72 = t64 * t66;
	t68 = cos(qJ(2));
	t71 = t64 * t68;
	t63 = cos(pkin(10));
	t57 = t60 * t71 + t63 * t66;
	t58 = -t60 * t72 + t63 * t68;
	t59 = sin(pkin(11));
	t62 = cos(pkin(11));
	t50 = t57 * t59 + t58 * t62;
	t65 = sin(qJ(5));
	t67 = cos(qJ(5));
	t44 = t50 * t67 - t65 * t73;
	t42 = 0.1e1 / t44 ^ 2;
	t43 = t50 * t65 + t67 * t73;
	t70 = t43 ^ 2 * t42 + 0.1e1;
	t49 = -t57 * t62 + t58 * t59;
	t69 = -t60 * t66 + t63 * t71;
	t56 = t60 * t68 + t63 * t72;
	t55 = (t59 * t68 - t62 * t66) * t61;
	t54 = (t59 * t66 + t62 * t68) * t61;
	t53 = 0.1e1 / t54 ^ 2;
	t46 = t56 * t59 + t69 * t62;
	t45 = -t56 * t62 + t69 * t59;
	t41 = atan2(-t46, t54);
	t39 = cos(t41);
	t38 = sin(t41);
	t37 = 0.1e1 / t70;
	t36 = -t38 * t46 + t39 * t54;
	t35 = 0.1e1 / t36 ^ 2;
	t33 = (-t45 / t54 + t55 * t46 * t53) / (t46 ^ 2 * t53 + 0.1e1);
	t1 = [0, t33, 0, 0, 0, 0; 0, (-t50 / t36 - (-t38 * t45 + t39 * t55 + (-t38 * t54 - t39 * t46) * t33) * t49 * t35) / (t49 ^ 2 * t35 + 0.1e1), 0, 0, 0, 0; 0, (t65 / t44 - t67 * t43 * t42) * t49 * t37, 0, 0, t70 * t37, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (632->36), mult. (1727->92), div. (65->9), fcn. (2425->15), ass. (0->56)
	t89 = sin(pkin(6));
	t97 = cos(qJ(5));
	t104 = t89 * t97;
	t92 = cos(pkin(6));
	t98 = cos(qJ(2));
	t102 = t92 * t98;
	t88 = sin(pkin(10));
	t91 = cos(pkin(10));
	t95 = sin(qJ(2));
	t82 = -t91 * t102 + t88 * t95;
	t103 = t92 * t95;
	t83 = t91 * t103 + t88 * t98;
	t87 = sin(pkin(11));
	t90 = cos(pkin(11));
	t70 = t82 * t87 + t83 * t90;
	t94 = sin(qJ(5));
	t64 = -t91 * t104 + t70 * t94;
	t81 = (-t87 * t98 + t90 * t95) * t89;
	t77 = t81 * t94 + t92 * t97;
	t63 = atan2(-t64, t77);
	t60 = sin(t63);
	t61 = cos(t63);
	t54 = -t60 * t64 + t61 * t77;
	t53 = 0.1e1 / t54 ^ 2;
	t84 = t88 * t102 + t91 * t95;
	t85 = -t88 * t103 + t91 * t98;
	t74 = t84 * t87 + t85 * t90;
	t67 = t88 * t104 + t74 * t94;
	t109 = t53 * t67;
	t100 = -t84 * t90 + t85 * t87;
	t105 = t89 * t94;
	t68 = -t88 * t105 + t74 * t97;
	t93 = sin(qJ(6));
	t96 = cos(qJ(6));
	t59 = t100 * t93 + t68 * t96;
	t57 = 0.1e1 / t59 ^ 2;
	t58 = -t100 * t96 + t68 * t93;
	t108 = t57 * t58;
	t76 = 0.1e1 / t77 ^ 2;
	t107 = t64 * t76;
	t106 = t100 * t97;
	t101 = t58 ^ 2 * t57 + 0.1e1;
	t99 = -t60 * t77 - t61 * t64;
	t80 = (t87 * t95 + t90 * t98) * t89;
	t78 = t81 * t97 - t92 * t94;
	t75 = 0.1e1 / t77;
	t69 = -t82 * t90 + t83 * t87;
	t66 = t91 * t105 + t70 * t97;
	t62 = 0.1e1 / (t64 ^ 2 * t76 + 0.1e1);
	t56 = 0.1e1 / t59;
	t55 = 0.1e1 / t101;
	t52 = 0.1e1 / t54;
	t51 = 0.1e1 / (t67 ^ 2 * t53 + 0.1e1);
	t50 = (t80 * t107 - t69 * t75) * t94 * t62;
	t49 = (t78 * t107 - t66 * t75) * t62;
	t1 = [0, t50, 0, 0, t49, 0; 0, (t100 * t94 * t52 - ((-t60 * t69 + t61 * t80) * t94 + t99 * t50) * t109) * t51, 0, 0, (t68 * t52 - (t99 * t49 - t60 * t66 + t61 * t78) * t109) * t51, 0; 0, ((t93 * t106 + t74 * t96) * t56 - (t96 * t106 - t74 * t93) * t108) * t55, 0, 0, (t96 * t108 - t56 * t93) * t67 * t55, t101 * t55;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end