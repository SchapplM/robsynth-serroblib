% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP10
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
% Datum: 2019-10-10 10:44
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP10_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP10_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:44:48
	% EndTime: 2019-10-10 10:44:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:44:48
	% EndTime: 2019-10-10 10:44:48
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
	% StartTime: 2019-10-10 10:44:48
	% EndTime: 2019-10-10 10:44:48
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
	% StartTime: 2019-10-10 10:44:48
	% EndTime: 2019-10-10 10:44:48
	% DurationCPUTime: 0.06s
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
	% StartTime: 2019-10-10 10:44:48
	% EndTime: 2019-10-10 10:44:48
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
	% StartTime: 2019-10-10 10:44:48
	% EndTime: 2019-10-10 10:44:48
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
	% StartTime: 2019-10-10 10:44:48
	% EndTime: 2019-10-10 10:44:48
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (122->30), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
	t174 = cos(pkin(6));
	t176 = sin(qJ(2));
	t180 = cos(qJ(1));
	t182 = t180 * t176;
	t177 = sin(qJ(1));
	t179 = cos(qJ(2));
	t184 = t177 * t179;
	t165 = t174 * t182 + t184;
	t172 = pkin(11) + qJ(4);
	t170 = sin(t172);
	t171 = cos(t172);
	t173 = sin(pkin(6));
	t187 = t173 * t180;
	t159 = -t165 * t171 + t170 * t187;
	t181 = t180 * t179;
	t185 = t177 * t176;
	t164 = -t174 * t181 + t185;
	t175 = sin(qJ(5));
	t178 = cos(qJ(5));
	t195 = t159 * t175 + t164 * t178;
	t194 = t159 * t178 - t164 * t175;
	t191 = t171 * t175;
	t190 = t171 * t178;
	t189 = t173 * t176;
	t188 = t173 * t177;
	t186 = t175 * t179;
	t183 = t178 * t179;
	t157 = -t165 * t170 - t171 * t187;
	t167 = -t174 * t185 + t181;
	t166 = t174 * t184 + t182;
	t163 = t174 * t170 + t171 * t189;
	t162 = -t170 * t189 + t174 * t171;
	t161 = t167 * t171 + t170 * t188;
	t160 = t167 * t170 - t171 * t188;
	t156 = t161 * t178 + t166 * t175;
	t155 = t161 * t175 - t166 * t178;
	t1 = [t194, -t166 * t190 + t167 * t175, 0, -t160 * t178, -t155, 0; t156, -t164 * t190 + t165 * t175, 0, t157 * t178, t195, 0; 0, (t171 * t183 + t175 * t176) * t173, 0, t162 * t178, -t163 * t175 - t173 * t183, 0; t157, -t166 * t170, 0, t161, 0, 0; t160, -t164 * t170, 0, -t159, 0, 0; 0, t173 * t179 * t170, 0, t163, 0, 0; t195, -t166 * t191 - t167 * t178, 0, -t160 * t175, t156, 0; t155, -t164 * t191 - t165 * t178, 0, t157 * t175, -t194, 0; 0, (t171 * t186 - t176 * t178) * t173, 0, t162 * t175, t163 * t178 - t173 * t186, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end