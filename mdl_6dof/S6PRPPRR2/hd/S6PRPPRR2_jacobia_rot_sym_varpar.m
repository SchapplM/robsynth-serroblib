% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPPRR2
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
%   Wie in S6PRPPRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:26
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPPRR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->0), mult. (78->0), div. (6->0), fcn. (108->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (161->15), mult. (461->36), div. (29->11), fcn. (658->11), ass. (0->26)
	t49 = cos(pkin(6));
	t44 = sin(pkin(11));
	t47 = cos(pkin(11));
	t50 = sin(qJ(2));
	t51 = cos(qJ(2));
	t53 = t51 * t44 + t50 * t47;
	t41 = t53 * t49;
	t42 = t50 * t44 - t51 * t47;
	t45 = sin(pkin(10));
	t48 = cos(pkin(10));
	t36 = -t45 * t41 - t48 * t42;
	t52 = t42 * t49;
	t46 = sin(pkin(6));
	t40 = t53 * t46;
	t39 = t42 * t46;
	t38 = 0.1e1 / t39 ^ 2;
	t34 = t45 * t52 - t48 * t53;
	t33 = -t45 * t53 - t48 * t52;
	t32 = -t48 * t41 + t45 * t42;
	t31 = atan2(t33, t39);
	t29 = cos(t31);
	t28 = sin(t31);
	t27 = t28 * t33 + t29 * t39;
	t26 = 0.1e1 / t27 ^ 2;
	t24 = (t32 / t39 - t40 * t33 * t38) / (t33 ^ 2 * t38 + 0.1e1);
	t1 = [0, t24, 0, 0, 0, 0; 0, (t36 / t27 + (t28 * t32 + t29 * t40 + (-t28 * t39 + t29 * t33) * t24) * t34 * t26) / (t34 ^ 2 * t26 + 0.1e1), 0, 0, 0, 0; 0, t34 / t45 / t46 / (0.1e1 + t36 ^ 2 / t45 ^ 2 / t46 ^ 2), 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (222->19), mult. (626->46), div. (35->9), fcn. (880->13), ass. (0->34)
	t59 = sin(pkin(10));
	t60 = sin(pkin(6));
	t71 = t59 * t60;
	t62 = cos(pkin(10));
	t58 = sin(pkin(11));
	t61 = cos(pkin(11));
	t65 = sin(qJ(2));
	t67 = cos(qJ(2));
	t55 = t58 * t65 - t67 * t61;
	t63 = cos(pkin(6));
	t68 = t55 * t63;
	t69 = t58 * t67 + t65 * t61;
	t48 = t59 * t68 - t62 * t69;
	t64 = sin(qJ(5));
	t66 = cos(qJ(5));
	t43 = -t48 * t64 + t66 * t71;
	t41 = 0.1e1 / t43 ^ 2;
	t42 = t48 * t66 + t64 * t71;
	t70 = t41 * t42 ^ 2 + 0.1e1;
	t54 = t69 * t63;
	t49 = -t54 * t59 - t55 * t62;
	t53 = t69 * t60;
	t52 = t55 * t60;
	t51 = 0.1e1 / t53 ^ 2;
	t45 = t54 * t62 - t55 * t59;
	t44 = -t59 * t69 - t62 * t68;
	t40 = atan2(-t45, t53);
	t38 = cos(t40);
	t37 = sin(t40);
	t36 = 0.1e1 / t70;
	t35 = -t37 * t45 + t38 * t53;
	t34 = 0.1e1 / t35 ^ 2;
	t32 = (-t44 / t53 - t52 * t45 * t51) / (t45 ^ 2 * t51 + 0.1e1);
	t1 = [0, t32, 0, 0, 0, 0; 0, (t48 / t35 - (-t37 * t44 - t38 * t52 + (-t37 * t53 - t38 * t45) * t32) * t49 * t34) / (t34 * t49 ^ 2 + 0.1e1), 0, 0, 0, 0; 0, (-t66 / t43 - t64 * t42 * t41) * t49 * t36, 0, 0, t70 * t36, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (632->31), mult. (1727->84), div. (65->9), fcn. (2425->15), ass. (0->54)
	t89 = sin(pkin(6));
	t94 = sin(qJ(5));
	t104 = t89 * t94;
	t87 = sin(pkin(11));
	t90 = cos(pkin(11));
	t95 = sin(qJ(2));
	t98 = cos(qJ(2));
	t85 = t95 * t87 - t98 * t90;
	t92 = cos(pkin(6));
	t83 = t85 * t92;
	t88 = sin(pkin(10));
	t91 = cos(pkin(10));
	t99 = t98 * t87 + t95 * t90;
	t73 = t91 * t83 + t88 * t99;
	t97 = cos(qJ(5));
	t69 = t91 * t104 + t73 * t97;
	t81 = t85 * t89;
	t79 = -t81 * t97 + t92 * t94;
	t66 = atan2(t69, t79);
	t63 = sin(t66);
	t64 = cos(t66);
	t57 = t63 * t69 + t64 * t79;
	t56 = 0.1e1 / t57 ^ 2;
	t75 = t88 * t83 - t91 * t99;
	t67 = t88 * t104 + t75 * t97;
	t108 = t56 * t67;
	t84 = t99 * t92;
	t100 = -t88 * t84 - t91 * t85;
	t103 = t89 * t97;
	t68 = t88 * t103 - t75 * t94;
	t93 = sin(qJ(6));
	t96 = cos(qJ(6));
	t62 = t100 * t93 + t68 * t96;
	t60 = 0.1e1 / t62 ^ 2;
	t61 = -t100 * t96 + t68 * t93;
	t107 = t60 * t61;
	t78 = 0.1e1 / t79 ^ 2;
	t106 = t69 * t78;
	t105 = t100 * t94;
	t102 = t61 ^ 2 * t60 + 0.1e1;
	t101 = -t63 * t79 + t64 * t69;
	t82 = t99 * t89;
	t80 = t81 * t94 + t92 * t97;
	t77 = 0.1e1 / t79;
	t72 = t91 * t84 - t88 * t85;
	t70 = t91 * t103 - t73 * t94;
	t65 = 0.1e1 / (t69 ^ 2 * t78 + 0.1e1);
	t59 = 0.1e1 / t62;
	t58 = 0.1e1 / t102;
	t55 = 0.1e1 / t57;
	t54 = 0.1e1 / (t67 ^ 2 * t56 + 0.1e1);
	t53 = (t82 * t106 + t72 * t77) * t97 * t65;
	t52 = (-t80 * t106 + t70 * t77) * t65;
	t1 = [0, t53, 0, 0, t52, 0; 0, (-t100 * t97 * t55 - ((t63 * t72 - t64 * t82) * t97 + t101 * t53) * t108) * t54, 0, 0, (t68 * t55 - (t101 * t52 + t63 * t70 + t64 * t80) * t108) * t54, 0; 0, ((t93 * t105 - t75 * t96) * t59 - (t96 * t105 + t75 * t93) * t107) * t58, 0, 0, (t96 * t107 - t59 * t93) * t67 * t58, t102 * t58;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end