% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRRRR15_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRRRR15_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR15_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
Ja_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (35->13), mult. (116->34), div. (25->9), fcn. (175->9), ass. (0->25)
	t30 = cos(pkin(5));
	t29 = sin(pkin(5));
	t34 = cos(qJ(1));
	t38 = t34 * t29;
	t26 = atan2(t38, t30);
	t23 = sin(t26);
	t24 = cos(t26);
	t18 = t23 * t38 + t24 * t30;
	t32 = sin(qJ(1));
	t42 = 0.1e1 / t18 ^ 2 * t32 ^ 2;
	t27 = t29 ^ 2;
	t25 = 0.1e1 / (0.1e1 + t34 ^ 2 * t27 / t30 ^ 2);
	t41 = t25 / t30;
	t31 = sin(qJ(2));
	t40 = t32 * t31;
	t33 = cos(qJ(2));
	t39 = t32 * t33;
	t37 = t34 * t31;
	t36 = t34 * t33;
	t22 = -t30 * t40 + t36;
	t20 = 0.1e1 / t22 ^ 2;
	t21 = t30 * t39 + t37;
	t35 = t21 ^ 2 * t20 + 0.1e1;
	t19 = 0.1e1 / t35;
	t1 = [-t32 * t29 * t41, 0, 0, 0, 0; (0.1e1 / t18 * t38 - (-t24 * t27 * t34 * t41 + (t25 - 0.1e1) * t29 * t23) * t29 * t42) / (t27 * t42 + 0.1e1), 0, 0, 0, 0; ((t30 * t36 - t40) / t22 - (-t30 * t37 - t39) * t21 * t20) * t19, t35 * t19, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (211->19), mult. (155->39), div. (32->9), fcn. (230->13), ass. (0->28)
	t58 = qJ(2) + qJ(3);
	t52 = pkin(5) + t58;
	t53 = pkin(5) - t58;
	t68 = -sin(t52) / 0.2e1 + sin(t53) / 0.2e1;
	t67 = cos(t53) / 0.2e1 + cos(t52) / 0.2e1;
	t60 = cos(pkin(5));
	t59 = sin(pkin(5));
	t62 = cos(qJ(1));
	t63 = t62 * t59;
	t50 = atan2(t63, t60);
	t47 = sin(t50);
	t48 = cos(t50);
	t40 = t47 * t63 + t48 * t60;
	t61 = sin(qJ(1));
	t66 = 0.1e1 / t40 ^ 2 * t61 ^ 2;
	t55 = cos(t58);
	t44 = t62 * t55 + t61 * t68;
	t42 = 0.1e1 / t44 ^ 2;
	t54 = sin(t58);
	t43 = t62 * t54 + t61 * t67;
	t65 = t43 ^ 2 * t42;
	t56 = t59 ^ 2;
	t49 = 0.1e1 / (0.1e1 + t62 ^ 2 * t56 / t60 ^ 2);
	t64 = t49 / t60;
	t41 = 0.1e1 / t44;
	t37 = 0.1e1 / (0.1e1 + t65);
	t36 = (t44 * t41 + t65) * t37;
	t1 = [-t61 * t59 * t64, 0, 0, 0, 0; (0.1e1 / t40 * t63 - (-t48 * t56 * t62 * t64 + (t49 - 0.1e1) * t59 * t47) * t59 * t66) / (t56 * t66 + 0.1e1), 0, 0, 0, 0; ((-t61 * t54 + t62 * t67) * t41 - (-t61 * t55 + t62 * t68) * t43 * t42) * t37, t36, t36, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (385->19), mult. (224->39), div. (38->9), fcn. (275->13), ass. (0->28)
	t71 = cos(pkin(5));
	t70 = sin(pkin(5));
	t73 = cos(qJ(1));
	t74 = t73 * t70;
	t61 = atan2(t74, t71);
	t58 = sin(t61);
	t59 = cos(t61);
	t51 = t58 * t74 + t59 * t71;
	t72 = sin(qJ(1));
	t77 = 0.1e1 / t51 ^ 2 * t72 ^ 2;
	t67 = qJ(2) + qJ(3) + qJ(4);
	t63 = pkin(5) + t67;
	t64 = pkin(5) - t67;
	t56 = sin(t63) / 0.2e1 - sin(t64) / 0.2e1;
	t66 = cos(t67);
	t55 = -t72 * t56 + t73 * t66;
	t53 = 0.1e1 / t55 ^ 2;
	t57 = cos(t64) / 0.2e1 + cos(t63) / 0.2e1;
	t65 = sin(t67);
	t54 = t72 * t57 + t73 * t65;
	t76 = t54 ^ 2 * t53;
	t68 = t70 ^ 2;
	t60 = 0.1e1 / (0.1e1 + t73 ^ 2 * t68 / t71 ^ 2);
	t75 = t60 / t71;
	t52 = 0.1e1 / t55;
	t48 = 0.1e1 / (0.1e1 + t76);
	t47 = (t55 * t52 + t76) * t48;
	t1 = [-t72 * t70 * t75, 0, 0, 0, 0; (0.1e1 / t51 * t74 - (-t59 * t68 * t73 * t75 + (t60 - 0.1e1) * t70 * t58) * t70 * t77) / (t68 * t77 + 0.1e1), 0, 0, 0, 0; ((t73 * t57 - t72 * t65) * t52 - (-t73 * t56 - t72 * t66) * t54 * t53) * t48, t47, t47, t47, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (1542->42), mult. (1953->93), div. (116->9), fcn. (2696->13), ass. (0->57)
	t125 = cos(pkin(6));
	t126 = cos(pkin(5));
	t128 = sin(qJ(1));
	t129 = cos(qJ(5));
	t140 = t129 * t128;
	t133 = t126 * t140;
	t127 = sin(qJ(5));
	t130 = cos(qJ(1));
	t139 = t130 * t127;
	t114 = -t125 * t133 - t139;
	t141 = t128 * t127;
	t134 = t126 * t141;
	t138 = t130 * t129;
	t116 = -t125 * t138 + t134;
	t122 = qJ(2) + qJ(3) + qJ(4);
	t120 = sin(t122);
	t121 = cos(t122);
	t123 = sin(pkin(6));
	t124 = sin(pkin(5));
	t143 = t124 * t128;
	t136 = t123 * t143;
	t101 = t114 * t121 + t116 * t120 + t129 * t136;
	t113 = t125 * t134 - t138;
	t115 = t125 * t139 + t133;
	t100 = -t113 * t121 - t115 * t120 + t127 * t136;
	t99 = 0.1e1 / t100 ^ 2;
	t149 = t101 * t99;
	t148 = t101 ^ 2 * t99;
	t142 = t128 * t121;
	t145 = t120 * t130;
	t107 = t125 * t143 + (t126 * t142 + t145) * t123;
	t146 = t120 * t128;
	t108 = -t123 * t146 + (t121 * t123 * t126 + t124 * t125) * t130;
	t144 = t123 * t124;
	t112 = -t121 * t144 + t126 * t125;
	t106 = atan2(t108, t112);
	t103 = sin(t106);
	t104 = cos(t106);
	t96 = t103 * t108 + t104 * t112;
	t95 = 0.1e1 / t96 ^ 2;
	t147 = t107 ^ 2 * t95;
	t137 = t120 * t144;
	t135 = t126 * t139;
	t132 = t130 * t144;
	t131 = t126 * t138;
	t111 = 0.1e1 / t112 ^ 2;
	t110 = 0.1e1 / t112;
	t109 = (t126 * t145 + t142) * t123;
	t105 = 0.1e1 / (t108 ^ 2 * t111 + 0.1e1);
	t98 = 0.1e1 / t100;
	t97 = 0.1e1 / (0.1e1 + t148);
	t94 = 0.1e1 / t96;
	t93 = 0.1e1 / (0.1e1 + t147);
	t92 = (-t108 * t111 * t137 - t109 * t110) * t105;
	t91 = ((t114 * t120 - t116 * t121) * t98 + (t113 * t120 - t115 * t121) * t149) * t97;
	t90 = ((t121 * t130 - t126 * t146) * t123 * t94 - ((t108 * t92 + t137) * t104 + (-t112 * t92 - t109) * t103) * t107 * t95) * t93;
	t1 = [-t107 * t110 * t105, t92, t92, t92, 0; (t108 * t94 - (-t103 + (-t104 * t108 * t110 + t103) * t105) * t147) * t93, t90, t90, t90, 0; ((-(-t125 * t131 + t141) * t121 - (t125 * t140 + t135) * t120 - t129 * t132) * t98 + ((-t125 * t135 - t140) * t121 + (t125 * t141 - t131) * t120 + t127 * t132) * t149) * t97, t91, t91, t91, (t100 * t98 + t148) * t97;];
	Ja_rot = t1;
end