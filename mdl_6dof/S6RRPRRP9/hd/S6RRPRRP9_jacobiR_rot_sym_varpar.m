% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:42
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP9_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP9_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:42:53
	% EndTime: 2019-10-10 10:42:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:42:53
	% EndTime: 2019-10-10 10:42:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:42:53
	% EndTime: 2019-10-10 10:42:53
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t50 = sin(qJ(2));
	t51 = sin(qJ(1));
	t57 = t51 * t50;
	t52 = cos(qJ(2));
	t56 = t51 * t52;
	t53 = cos(qJ(1));
	t55 = t53 * t50;
	t54 = t53 * t52;
	t49 = cos(pkin(6));
	t48 = sin(pkin(6));
	t47 = -t49 * t57 + t54;
	t46 = -t49 * t56 - t55;
	t45 = -t49 * t55 - t56;
	t44 = -t49 * t54 + t57;
	t1 = [t45, t46, 0, 0, 0, 0; t47, -t44, 0, 0, 0, 0; 0, t48 * t52, 0, 0, 0, 0; t44, -t47, 0, 0, 0, 0; t46, t45, 0, 0, 0, 0; 0, -t48 * t50, 0, 0, 0, 0; t53 * t48, 0, 0, 0, 0, 0; t51 * t48, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:42:53
	% EndTime: 2019-10-10 10:42:53
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (19->11), mult. (57->26), div. (0->0), fcn. (88->8), ass. (0->20)
	t66 = sin(pkin(6));
	t70 = sin(qJ(1));
	t79 = t66 * t70;
	t71 = cos(qJ(2));
	t78 = t66 * t71;
	t72 = cos(qJ(1));
	t77 = t66 * t72;
	t69 = sin(qJ(2));
	t76 = t70 * t69;
	t75 = t70 * t71;
	t74 = t72 * t69;
	t73 = t72 * t71;
	t68 = cos(pkin(6));
	t67 = cos(pkin(11));
	t65 = sin(pkin(11));
	t64 = -t68 * t76 + t73;
	t63 = t68 * t75 + t74;
	t62 = t68 * t74 + t75;
	t61 = t68 * t73 - t76;
	t1 = [-t62 * t67 + t65 * t77, -t63 * t67, 0, 0, 0, 0; t64 * t67 + t65 * t79, t61 * t67, 0, 0, 0, 0; 0, t67 * t78, 0, 0, 0, 0; t62 * t65 + t67 * t77, t63 * t65, 0, 0, 0, 0; -t64 * t65 + t67 * t79, -t61 * t65, 0, 0, 0, 0; 0, -t65 * t78, 0, 0, 0, 0; t61, t64, 0, 0, 0, 0; t63, t62, 0, 0, 0, 0; 0, t66 * t69, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:42:53
	% EndTime: 2019-10-10 10:42:53
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (55->16), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->26)
	t90 = sin(pkin(6));
	t92 = sin(qJ(2));
	t105 = t90 * t92;
	t93 = sin(qJ(1));
	t104 = t90 * t93;
	t94 = cos(qJ(2));
	t103 = t90 * t94;
	t95 = cos(qJ(1));
	t102 = t90 * t95;
	t101 = t93 * t92;
	t100 = t93 * t94;
	t99 = t95 * t92;
	t98 = t95 * t94;
	t91 = cos(pkin(6));
	t83 = t91 * t99 + t100;
	t89 = pkin(11) + qJ(4);
	t87 = sin(t89);
	t88 = cos(t89);
	t97 = t87 * t102 - t83 * t88;
	t96 = t88 * t102 + t83 * t87;
	t85 = -t91 * t101 + t98;
	t84 = t91 * t100 + t99;
	t82 = t91 * t98 - t101;
	t81 = t87 * t104 + t85 * t88;
	t80 = t88 * t104 - t85 * t87;
	t1 = [t97, -t84 * t88, 0, t80, 0, 0; t81, t82 * t88, 0, -t96, 0, 0; 0, t88 * t103, 0, -t87 * t105 + t91 * t88, 0, 0; t96, t84 * t87, 0, -t81, 0, 0; t80, -t82 * t87, 0, t97, 0, 0; 0, -t87 * t103, 0, -t88 * t105 - t91 * t87, 0, 0; t82, t85, 0, 0, 0, 0; t84, t83, 0, 0, 0, 0; 0, t105, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:42:54
	% EndTime: 2019-10-10 10:42:54
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (125->31), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
	t142 = cos(pkin(6));
	t144 = sin(qJ(2));
	t148 = cos(qJ(1));
	t150 = t148 * t144;
	t145 = sin(qJ(1));
	t147 = cos(qJ(2));
	t152 = t145 * t147;
	t133 = t142 * t150 + t152;
	t140 = pkin(11) + qJ(4);
	t138 = sin(t140);
	t139 = cos(t140);
	t141 = sin(pkin(6));
	t155 = t141 * t148;
	t127 = -t133 * t139 + t138 * t155;
	t149 = t148 * t147;
	t153 = t145 * t144;
	t132 = -t142 * t149 + t153;
	t143 = sin(qJ(5));
	t146 = cos(qJ(5));
	t163 = t127 * t143 + t132 * t146;
	t162 = t127 * t146 - t132 * t143;
	t159 = t139 * t143;
	t158 = t139 * t146;
	t157 = t141 * t144;
	t156 = t141 * t145;
	t154 = t143 * t147;
	t151 = t146 * t147;
	t125 = -t133 * t138 - t139 * t155;
	t135 = -t142 * t153 + t149;
	t134 = t142 * t152 + t150;
	t131 = t142 * t138 + t139 * t157;
	t130 = -t138 * t157 + t142 * t139;
	t129 = t135 * t139 + t138 * t156;
	t128 = t135 * t138 - t139 * t156;
	t124 = t129 * t146 + t134 * t143;
	t123 = -t129 * t143 + t134 * t146;
	t1 = [t162, -t134 * t158 + t135 * t143, 0, -t128 * t146, t123, 0; t124, -t132 * t158 + t133 * t143, 0, t125 * t146, t163, 0; 0, (t139 * t151 + t143 * t144) * t141, 0, t130 * t146, -t131 * t143 - t141 * t151, 0; -t163, t134 * t159 + t135 * t146, 0, t128 * t143, -t124, 0; t123, t132 * t159 + t133 * t146, 0, -t125 * t143, t162, 0; 0, (-t139 * t154 + t144 * t146) * t141, 0, -t130 * t143, -t131 * t146 + t141 * t154, 0; t125, -t134 * t138, 0, t129, 0, 0; t128, -t132 * t138, 0, -t127, 0, 0; 0, t141 * t147 * t138, 0, t131, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:42:54
	% EndTime: 2019-10-10 10:42:54
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (125->31), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
	t144 = cos(pkin(6));
	t146 = sin(qJ(2));
	t150 = cos(qJ(1));
	t152 = t150 * t146;
	t147 = sin(qJ(1));
	t149 = cos(qJ(2));
	t154 = t147 * t149;
	t135 = t144 * t152 + t154;
	t142 = pkin(11) + qJ(4);
	t140 = sin(t142);
	t141 = cos(t142);
	t143 = sin(pkin(6));
	t157 = t143 * t150;
	t129 = -t135 * t141 + t140 * t157;
	t151 = t150 * t149;
	t155 = t147 * t146;
	t134 = -t144 * t151 + t155;
	t145 = sin(qJ(5));
	t148 = cos(qJ(5));
	t165 = t129 * t145 + t134 * t148;
	t164 = t129 * t148 - t134 * t145;
	t161 = t141 * t145;
	t160 = t141 * t148;
	t159 = t143 * t146;
	t158 = t143 * t147;
	t156 = t145 * t149;
	t153 = t148 * t149;
	t127 = -t135 * t140 - t141 * t157;
	t137 = -t144 * t155 + t151;
	t136 = t144 * t154 + t152;
	t133 = t144 * t140 + t141 * t159;
	t132 = -t140 * t159 + t144 * t141;
	t131 = t137 * t141 + t140 * t158;
	t130 = t137 * t140 - t141 * t158;
	t126 = t131 * t148 + t136 * t145;
	t125 = -t131 * t145 + t136 * t148;
	t1 = [t164, -t136 * t160 + t137 * t145, 0, -t130 * t148, t125, 0; t126, -t134 * t160 + t135 * t145, 0, t127 * t148, t165, 0; 0, (t141 * t153 + t145 * t146) * t143, 0, t132 * t148, -t133 * t145 - t143 * t153, 0; -t165, t136 * t161 + t137 * t148, 0, t130 * t145, -t126, 0; t125, t134 * t161 + t135 * t148, 0, -t127 * t145, t164, 0; 0, (-t141 * t156 + t146 * t148) * t143, 0, -t132 * t145, -t133 * t148 + t143 * t156, 0; t127, -t136 * t140, 0, t131, 0, 0; t130, -t134 * t140, 0, -t129, 0, 0; 0, t143 * t149 * t140, 0, t133, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end