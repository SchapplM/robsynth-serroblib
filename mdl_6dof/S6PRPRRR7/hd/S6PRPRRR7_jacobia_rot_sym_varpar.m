% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR7
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
%   Wie in S6PRPRRR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:05
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobia_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (148->22), mult. (435->55), div. (30->9), fcn. (617->13), ass. (0->34)
	t43 = sin(pkin(13));
	t47 = cos(pkin(13));
	t51 = cos(qJ(2));
	t49 = cos(pkin(6));
	t50 = sin(qJ(2));
	t54 = t49 * t50;
	t40 = -t43 * t54 + t47 * t51;
	t48 = cos(pkin(7));
	t58 = t40 * t48;
	t44 = sin(pkin(7));
	t45 = sin(pkin(6));
	t57 = t44 * t45;
	t56 = t45 * t48;
	t55 = t45 * t50;
	t53 = t49 * t51;
	t39 = -t43 * t53 - t47 * t50;
	t52 = t39 * t48 + t43 * t57;
	t46 = cos(pkin(14));
	t42 = sin(pkin(14));
	t38 = -t43 * t51 - t47 * t54;
	t37 = t49 * t48 - t51 * t57;
	t36 = 0.1e1 / t37 ^ 2;
	t35 = -t39 * t44 + t43 * t56;
	t34 = (-t43 * t50 + t47 * t53) * t44 + t47 * t56;
	t33 = atan2(t34, t37);
	t31 = cos(t33);
	t30 = sin(t33);
	t29 = t40 * t46 + t52 * t42;
	t28 = t40 * t42 - t52 * t46;
	t27 = 0.1e1 / t29 ^ 2;
	t25 = t30 * t34 + t31 * t37;
	t24 = 0.1e1 / t25 ^ 2;
	t22 = (t38 / t37 - t34 * t36 * t55) * t44 / (t34 ^ 2 * t36 + 0.1e1);
	t1 = [0, t22, 0, 0, 0, 0; 0, (t40 * t44 / t25 - ((t30 * t38 + t31 * t55) * t44 + (-t30 * t37 + t31 * t34) * t22) * t35 * t24) / (t35 ^ 2 * t24 + 0.1e1), 0, 0, 0, 0; 0, ((t39 * t42 + t46 * t58) / t29 - (t39 * t46 - t42 * t58) * t28 * t27) / (t28 ^ 2 * t27 + 0.1e1), 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (563->40), mult. (1633->92), div. (35->9), fcn. (2228->17), ass. (0->52)
	t100 = cos(qJ(2));
	t93 = cos(pkin(13));
	t105 = t93 * t100;
	t96 = cos(pkin(6));
	t98 = sin(qJ(2));
	t107 = t96 * t98;
	t88 = sin(pkin(13));
	t85 = -t88 * t107 + t105;
	t87 = sin(pkin(14));
	t112 = t85 * t87;
	t90 = sin(pkin(7));
	t91 = sin(pkin(6));
	t111 = t90 * t91;
	t94 = cos(pkin(8));
	t110 = t90 * t94;
	t95 = cos(pkin(7));
	t109 = t91 * t95;
	t92 = cos(pkin(14));
	t108 = t92 * t95;
	t106 = t88 * t100;
	t84 = -t96 * t106 - t93 * t98;
	t101 = t88 * t111 + t84 * t95;
	t75 = t101 * t92 - t112;
	t81 = t88 * t109 - t84 * t90;
	t89 = sin(pkin(8));
	t103 = t75 * t94 + t81 * t89;
	t76 = t101 * t87 + t85 * t92;
	t97 = sin(qJ(4));
	t99 = cos(qJ(4));
	t65 = t103 * t97 + t76 * t99;
	t63 = 0.1e1 / t65 ^ 2;
	t64 = -t103 * t99 + t76 * t97;
	t104 = t64 ^ 2 * t63 + 0.1e1;
	t77 = -t85 * t108 - t84 * t87;
	t102 = t85 * t89 * t90 + t77 * t94;
	t83 = t93 * t107 + t106;
	t82 = t96 * t105 - t88 * t98;
	t79 = (-(-t100 * t87 - t98 * t108) * t89 + t98 * t110) * t91;
	t78 = -t95 * t112 + t84 * t92;
	t74 = -(t96 * t90 * t92 + (t100 * t108 - t87 * t98) * t91) * t89 + (-t100 * t111 + t96 * t95) * t94;
	t73 = 0.1e1 / t74 ^ 2;
	t72 = (-t83 * t108 - t82 * t87) * t89 - t83 * t110;
	t71 = -t75 * t89 + t81 * t94;
	t70 = (-t83 * t87 + (-t93 * t111 + t82 * t95) * t92) * t89 - (-t93 * t109 - t82 * t90) * t94;
	t69 = atan2(t70, t74);
	t67 = cos(t69);
	t66 = sin(t69);
	t62 = 0.1e1 / t104;
	t61 = t66 * t70 + t67 * t74;
	t60 = 0.1e1 / t61 ^ 2;
	t58 = (t72 / t74 - t79 * t70 * t73) / (t70 ^ 2 * t73 + 0.1e1);
	t1 = [0, t58, 0, 0, 0, 0; 0, ((t85 * t110 - t77 * t89) / t61 - (t66 * t72 + t67 * t79 + (-t66 * t74 + t67 * t70) * t58) * t71 * t60) / (t71 ^ 2 * t60 + 0.1e1), 0, 0, 0, 0; 0, ((-t102 * t99 + t78 * t97) / t65 - (t102 * t97 + t78 * t99) * t64 * t63) * t62, 0, t104 * t62, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:03
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (1772->61), mult. (5223->143), div. (65->9), fcn. (7061->19), ass. (0->78)
	t127 = sin(pkin(13));
	t132 = cos(pkin(13));
	t141 = cos(qJ(2));
	t135 = cos(pkin(6));
	t138 = sin(qJ(2));
	t154 = t135 * t138;
	t122 = t127 * t141 + t132 * t154;
	t126 = sin(pkin(14));
	t128 = sin(pkin(8));
	t129 = sin(pkin(7));
	t131 = cos(pkin(14));
	t133 = cos(pkin(8));
	t134 = cos(pkin(7));
	t153 = t135 * t141;
	t145 = -t127 * t138 + t132 * t153;
	t130 = sin(pkin(6));
	t158 = t130 * t129;
	t144 = -t132 * t158 + t145 * t134;
	t157 = t130 * t134;
	t168 = (-t122 * t126 + t144 * t131) * t133 + (-t145 * t129 - t132 * t157) * t128;
	t124 = -t127 * t154 + t132 * t141;
	t123 = -t127 * t153 - t132 * t138;
	t146 = t123 * t134 + t127 * t158;
	t112 = -t124 * t126 + t146 * t131;
	t120 = -t123 * t129 + t127 * t157;
	t167 = t112 * t133 + t120 * t128;
	t151 = t138 * t131;
	t155 = t134 * t141;
	t159 = t129 * t135;
	t119 = t130 * t151 + (t130 * t155 + t159) * t126;
	t137 = sin(qJ(4));
	t140 = cos(qJ(4));
	t152 = t138 * t126;
	t148 = (t131 * t159 + (t131 * t155 - t152) * t130) * t133 + (t135 * t134 - t141 * t158) * t128;
	t104 = t119 * t137 - t148 * t140;
	t111 = t122 * t131 + t144 * t126;
	t95 = t111 * t137 - t168 * t140;
	t94 = atan2(-t95, t104);
	t91 = sin(t94);
	t92 = cos(t94);
	t85 = t92 * t104 - t91 * t95;
	t84 = 0.1e1 / t85 ^ 2;
	t113 = t124 * t131 + t146 * t126;
	t98 = t113 * t137 - t167 * t140;
	t166 = t84 * t98;
	t106 = -t112 * t128 + t120 * t133;
	t136 = sin(qJ(5));
	t139 = cos(qJ(5));
	t99 = t113 * t140 + t167 * t137;
	t90 = t106 * t136 + t99 * t139;
	t88 = 0.1e1 / t90 ^ 2;
	t89 = -t106 * t139 + t99 * t136;
	t165 = t88 * t89;
	t103 = 0.1e1 / t104 ^ 2;
	t164 = t103 * t95;
	t161 = t126 * t134;
	t160 = t129 * t128;
	t156 = t134 * t131;
	t150 = t89 ^ 2 * t88 + 0.1e1;
	t149 = -t104 * t91 - t92 * t95;
	t114 = -t123 * t126 - t124 * t156;
	t147 = t114 * t133 + t124 * t160;
	t115 = t123 * t131 - t124 * t161;
	t108 = ((t131 * t141 - t134 * t152) * t137 + (-(-t126 * t141 - t134 * t151) * t133 - t138 * t160) * t140) * t130;
	t107 = t124 * t129 * t133 - t114 * t128;
	t105 = t119 * t140 + t148 * t137;
	t102 = 0.1e1 / t104;
	t101 = t115 * t140 + t147 * t137;
	t100 = (-t122 * t161 + t145 * t131) * t137 + (-(-t122 * t156 - t145 * t126) * t133 - t122 * t160) * t140;
	t97 = t111 * t140 + t168 * t137;
	t93 = 0.1e1 / (t95 ^ 2 * t103 + 0.1e1);
	t87 = 0.1e1 / t90;
	t86 = 0.1e1 / t150;
	t83 = 0.1e1 / t85;
	t82 = 0.1e1 / (t98 ^ 2 * t84 + 0.1e1);
	t81 = (-t100 * t102 + t108 * t164) * t93;
	t80 = (-t102 * t97 + t105 * t164) * t93;
	t1 = [0, t81, 0, t80, 0, 0; 0, ((t115 * t137 - t147 * t140) * t83 - (-t91 * t100 + t92 * t108 + t149 * t81) * t166) * t82, 0, (t99 * t83 - (t92 * t105 + t149 * t80 - t91 * t97) * t166) * t82, 0, 0; 0, ((t101 * t136 - t107 * t139) * t87 - (t101 * t139 + t107 * t136) * t165) * t86, 0, (-t136 * t87 + t139 * t165) * t98 * t86, t150 * t86, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:04
	% DurationCPUTime: 1.12s
	% Computational Cost: add. (4256->82), mult. (12302->195), div. (95->9), fcn. (16658->21), ass. (0->98)
	t166 = sin(pkin(8));
	t171 = cos(pkin(8));
	t165 = sin(pkin(13));
	t170 = cos(pkin(13));
	t181 = cos(qJ(2));
	t173 = cos(pkin(6));
	t177 = sin(qJ(2));
	t198 = t173 * t177;
	t162 = -t165 * t198 + t170 * t181;
	t164 = sin(pkin(14));
	t169 = cos(pkin(14));
	t197 = t173 * t181;
	t161 = -t165 * t197 - t170 * t177;
	t172 = cos(pkin(7));
	t167 = sin(pkin(7));
	t168 = sin(pkin(6));
	t206 = t167 * t168;
	t189 = t161 * t172 + t165 * t206;
	t186 = t162 * t164 - t189 * t169;
	t202 = t168 * t172;
	t190 = -t161 * t167 + t165 * t202;
	t213 = -t190 * t166 + t186 * t171;
	t160 = t165 * t181 + t170 * t198;
	t159 = -t165 * t177 + t170 * t197;
	t191 = t159 * t172 - t170 * t206;
	t149 = t160 * t169 + t191 * t164;
	t176 = sin(qJ(4));
	t180 = cos(qJ(4));
	t187 = -t160 * t164 + t191 * t169;
	t192 = -t159 * t167 - t170 * t202;
	t183 = t192 * t166 + t187 * t171;
	t136 = t149 * t180 + t183 * t176;
	t175 = sin(qJ(5));
	t179 = cos(qJ(5));
	t184 = -t187 * t166 + t192 * t171;
	t124 = t136 * t175 - t184 * t179;
	t199 = t172 * t181;
	t204 = t167 * t173;
	t156 = t168 * t177 * t169 + (t168 * t199 + t204) * t164;
	t155 = t169 * t204 + (-t164 * t177 + t169 * t199) * t168;
	t158 = t173 * t172 - t181 * t206;
	t194 = t155 * t171 + t158 * t166;
	t145 = t156 * t180 + t194 * t176;
	t148 = -t155 * t166 + t158 * t171;
	t133 = t145 * t175 - t148 * t179;
	t123 = atan2(-t124, t133);
	t120 = sin(t123);
	t121 = cos(t123);
	t114 = -t120 * t124 + t121 * t133;
	t113 = 0.1e1 / t114 ^ 2;
	t150 = t162 * t169 + t189 * t164;
	t138 = t150 * t180 - t213 * t176;
	t182 = t186 * t166 + t190 * t171;
	t127 = t138 * t175 - t182 * t179;
	t212 = t113 * t127;
	t128 = t138 * t179 + t182 * t175;
	t137 = t150 * t176 + t213 * t180;
	t174 = sin(qJ(6));
	t178 = cos(qJ(6));
	t119 = t128 * t178 + t137 * t174;
	t117 = 0.1e1 / t119 ^ 2;
	t118 = t128 * t174 - t137 * t178;
	t211 = t117 * t118;
	t132 = 0.1e1 / t133 ^ 2;
	t210 = t124 * t132;
	t209 = t137 * t179;
	t208 = t164 * t172;
	t207 = t167 * t166;
	t205 = t167 * t171;
	t203 = t167 * t177;
	t201 = t169 * t172;
	t200 = t172 * t177;
	t196 = t118 ^ 2 * t117 + 0.1e1;
	t195 = -t120 * t133 - t121 * t124;
	t152 = -t161 * t164 - t162 * t201;
	t193 = t152 * t171 + t162 * t207;
	t153 = t161 * t169 - t162 * t208;
	t151 = -t159 * t164 - t160 * t201;
	t146 = -t152 * t166 + t162 * t205;
	t144 = -t156 * t176 + t194 * t180;
	t141 = ((t171 * t176 * t175 + t166 * t179) * (-t164 * t181 - t169 * t200) + ((-t164 * t200 + t169 * t181) * t180 + t166 * t176 * t203) * t175 - t171 * t179 * t203) * t168;
	t140 = t153 * t180 + t193 * t176;
	t139 = t153 * t176 - t193 * t180;
	t135 = -t149 * t176 + t183 * t180;
	t134 = t145 * t179 + t148 * t175;
	t131 = 0.1e1 / t133;
	t130 = t140 * t179 + t146 * t175;
	t129 = ((t159 * t169 - t160 * t208) * t180 + (t151 * t171 + t160 * t207) * t176) * t175 - (-t151 * t166 + t160 * t205) * t179;
	t126 = t136 * t179 + t184 * t175;
	t122 = 0.1e1 / (t124 ^ 2 * t132 + 0.1e1);
	t116 = 0.1e1 / t119;
	t115 = 0.1e1 / t196;
	t112 = 0.1e1 / t114;
	t111 = 0.1e1 / (t127 ^ 2 * t113 + 0.1e1);
	t110 = (-t131 * t135 + t144 * t210) * t175 * t122;
	t109 = (-t129 * t131 + t141 * t210) * t122;
	t108 = (-t126 * t131 + t134 * t210) * t122;
	t1 = [0, t109, 0, t110, t108, 0; 0, ((t140 * t175 - t146 * t179) * t112 - (t195 * t109 - t120 * t129 + t121 * t141) * t212) * t111, 0, (-t137 * t175 * t112 - ((-t120 * t135 + t121 * t144) * t175 + t195 * t110) * t212) * t111, (t128 * t112 - (t195 * t108 - t120 * t126 + t121 * t134) * t212) * t111, 0; 0, ((t130 * t174 - t139 * t178) * t116 - (t130 * t178 + t139 * t174) * t211) * t115, 0, ((-t138 * t178 - t174 * t209) * t116 - (t138 * t174 - t178 * t209) * t211) * t115, (-t174 * t116 + t178 * t211) * t127 * t115, t196 * t115;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end