% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:01
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(11));
	t20 = sin(pkin(6));
	t19 = sin(pkin(11));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t46 = cos(pkin(6));
	t47 = sin(qJ(2));
	t50 = t46 * t47;
	t48 = cos(qJ(2));
	t49 = t46 * t48;
	t45 = cos(pkin(11));
	t44 = sin(pkin(6));
	t43 = sin(pkin(11));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, t43 * t49 + t45 * t47, 0, 0, 0, 0; 0, t43 * t47 - t45 * t49, 0, 0, 0, 0; 0, -t44 * t48, 0, 0, 0, 0; 0, -t43 * t50 + t45 * t48, 0, 0, 0, 0; 0, t43 * t48 + t45 * t50, 0, 0, 0, 0; 0, t44 * t47, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (16->12), mult. (57->31), div. (0->0), fcn. (88->8), ass. (0->18)
	t66 = sin(pkin(6));
	t69 = sin(qJ(4));
	t77 = t66 * t69;
	t71 = cos(qJ(4));
	t76 = t66 * t71;
	t72 = cos(qJ(2));
	t75 = t66 * t72;
	t68 = cos(pkin(6));
	t70 = sin(qJ(2));
	t74 = t68 * t70;
	t73 = t68 * t72;
	t67 = cos(pkin(11));
	t65 = sin(pkin(11));
	t64 = -t65 * t74 + t67 * t72;
	t63 = t65 * t73 + t67 * t70;
	t62 = t65 * t72 + t67 * t74;
	t61 = t65 * t70 - t67 * t73;
	t1 = [0, t64 * t69, 0, t63 * t71 - t65 * t77, 0, 0; 0, t62 * t69, 0, t61 * t71 + t67 * t77, 0, 0; 0, t70 * t77, 0, -t68 * t69 - t71 * t75, 0, 0; 0, t64 * t71, 0, -t63 * t69 - t65 * t76, 0, 0; 0, t62 * t71, 0, -t61 * t69 + t67 * t76, 0, 0; 0, t70 * t76, 0, -t68 * t71 + t69 * t75, 0, 0; 0, -t63, 0, 0, 0, 0; 0, -t61, 0, 0, 0, 0; 0, t75, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:17
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (56->13), mult. (87->32), div. (0->0), fcn. (134->8), ass. (0->26)
	t92 = sin(pkin(11));
	t93 = sin(pkin(6));
	t103 = t92 * t93;
	t94 = cos(pkin(11));
	t102 = t93 * t94;
	t96 = sin(qJ(2));
	t101 = t93 * t96;
	t97 = cos(qJ(2));
	t100 = t93 * t97;
	t95 = cos(pkin(6));
	t99 = t95 * t96;
	t98 = t95 * t97;
	t91 = qJ(4) + qJ(5);
	t90 = cos(t91);
	t89 = sin(t91);
	t88 = -t92 * t99 + t94 * t97;
	t87 = t92 * t98 + t94 * t96;
	t86 = t92 * t97 + t94 * t99;
	t85 = t92 * t96 - t94 * t98;
	t84 = t89 * t100 - t95 * t90;
	t83 = -t90 * t100 - t95 * t89;
	t82 = t90 * t102 - t85 * t89;
	t81 = t89 * t102 + t85 * t90;
	t80 = -t90 * t103 - t87 * t89;
	t79 = -t89 * t103 + t87 * t90;
	t1 = [0, t88 * t89, 0, t79, t79, 0; 0, t86 * t89, 0, t81, t81, 0; 0, t89 * t101, 0, t83, t83, 0; 0, t88 * t90, 0, t80, t80, 0; 0, t86 * t90, 0, t82, t82, 0; 0, t90 * t101, 0, t84, t84, 0; 0, -t87, 0, 0, 0, 0; 0, -t85, 0, 0, 0, 0; 0, t100, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:17
	% EndTime: 2019-10-09 22:01:17
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (132->32), mult. (214->65), div. (0->0), fcn. (313->10), ass. (0->37)
	t140 = sin(pkin(11));
	t142 = cos(pkin(11));
	t145 = sin(qJ(2));
	t143 = cos(pkin(6));
	t147 = cos(qJ(2));
	t150 = t143 * t147;
	t134 = t140 * t150 + t142 * t145;
	t139 = qJ(4) + qJ(5);
	t137 = sin(t139);
	t138 = cos(t139);
	t141 = sin(pkin(6));
	t154 = t140 * t141;
	t125 = t134 * t138 - t137 * t154;
	t144 = sin(qJ(6));
	t159 = t125 * t144;
	t132 = t140 * t145 - t142 * t150;
	t153 = t141 * t142;
	t127 = t132 * t138 + t137 * t153;
	t158 = t127 * t144;
	t152 = t141 * t147;
	t130 = -t143 * t137 - t138 * t152;
	t157 = t130 * t144;
	t156 = t137 * t144;
	t146 = cos(qJ(6));
	t155 = t137 * t146;
	t151 = t143 * t145;
	t149 = t144 * t145;
	t148 = t145 * t146;
	t135 = -t140 * t151 + t142 * t147;
	t133 = t140 * t147 + t142 * t151;
	t131 = -t137 * t152 + t143 * t138;
	t129 = t130 * t146;
	t128 = t132 * t137 - t138 * t153;
	t126 = t134 * t137 + t138 * t154;
	t124 = t127 * t146;
	t123 = t125 * t146;
	t1 = [0, -t134 * t144 + t135 * t155, 0, t123, t123, -t126 * t144 + t135 * t146; 0, -t132 * t144 + t133 * t155, 0, t124, t124, -t128 * t144 + t133 * t146; 0, (t137 * t148 + t144 * t147) * t141, 0, t129, t129, -t131 * t144 + t141 * t148; 0, -t134 * t146 - t135 * t156, 0, -t159, -t159, -t126 * t146 - t135 * t144; 0, -t132 * t146 - t133 * t156, 0, -t158, -t158, -t128 * t146 - t133 * t144; 0, (-t137 * t149 + t146 * t147) * t141, 0, -t157, -t157, -t131 * t146 - t141 * t149; 0, -t135 * t138, 0, t126, t126, 0; 0, -t133 * t138, 0, t128, t128, 0; 0, -t141 * t145 * t138, 0, t131, t131, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end