% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR6
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
%   Wie in S6PRRRRR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:23
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRR6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobia_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:36
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (179->22), mult. (525->56), div. (35->9), fcn. (736->13), ass. (0->37)
	t45 = sin(pkin(14));
	t48 = cos(pkin(14));
	t54 = cos(qJ(2));
	t50 = cos(pkin(6));
	t52 = sin(qJ(2));
	t58 = t50 * t52;
	t43 = -t45 * t58 + t48 * t54;
	t51 = sin(qJ(3));
	t63 = t43 * t51;
	t53 = cos(qJ(3));
	t62 = t43 * t53;
	t46 = sin(pkin(7));
	t47 = sin(pkin(6));
	t61 = t46 * t47;
	t49 = cos(pkin(7));
	t60 = t47 * t49;
	t59 = t47 * t52;
	t57 = t50 * t54;
	t42 = -t45 * t57 - t48 * t52;
	t55 = t42 * t49 + t45 * t61;
	t32 = t51 * t55 + t62;
	t30 = 0.1e1 / t32 ^ 2;
	t31 = -t53 * t55 + t63;
	t56 = t31 ^ 2 * t30 + 0.1e1;
	t41 = -t45 * t54 - t48 * t58;
	t40 = t50 * t49 - t54 * t61;
	t39 = 0.1e1 / t40 ^ 2;
	t38 = -t42 * t46 + t45 * t60;
	t37 = (-t45 * t52 + t48 * t57) * t46 + t48 * t60;
	t36 = atan2(t37, t40);
	t34 = cos(t36);
	t33 = sin(t36);
	t29 = 0.1e1 / t56;
	t28 = t33 * t37 + t34 * t40;
	t27 = 0.1e1 / t28 ^ 2;
	t25 = (t41 / t40 - t37 * t39 * t59) * t46 / (t37 ^ 2 * t39 + 0.1e1);
	t1 = [0, t25, 0, 0, 0, 0; 0, (t43 * t46 / t28 - ((t33 * t41 + t34 * t59) * t46 + (-t33 * t40 + t34 * t37) * t25) * t38 * t27) / (t38 ^ 2 * t27 + 0.1e1), 0, 0, 0, 0; 0, ((t42 * t51 + t49 * t62) / t32 - (t42 * t53 - t49 * t63) * t31 * t30) * t29, t56 * t29, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:36
	% EndTime: 2019-10-09 23:23:36
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (1034->50), mult. (3001->123), div. (65->9), fcn. (4101->17), ass. (0->69)
	t103 = sin(pkin(8));
	t104 = sin(pkin(7));
	t106 = cos(pkin(14));
	t107 = cos(pkin(8));
	t111 = sin(qJ(3));
	t114 = cos(qJ(3));
	t108 = cos(pkin(7));
	t105 = sin(pkin(6));
	t130 = t104 * t105;
	t102 = sin(pkin(14));
	t112 = sin(qJ(2));
	t109 = cos(pkin(6));
	t115 = cos(qJ(2));
	t124 = t109 * t115;
	t97 = -t102 * t112 + t106 * t124;
	t116 = t106 * t130 - t108 * t97;
	t127 = t105 * t108;
	t125 = t109 * t112;
	t98 = t102 * t115 + t106 * t125;
	t82 = (-t98 * t111 - t116 * t114) * t103 - (-t97 * t104 - t106 * t127) * t107;
	t126 = t108 * t114;
	t128 = t104 * t109;
	t87 = -(t114 * t128 + (-t111 * t112 + t115 * t126) * t105) * t103 + (t109 * t108 - t115 * t130) * t107;
	t81 = atan2(t82, t87);
	t78 = sin(t81);
	t79 = cos(t81);
	t72 = t78 * t82 + t79 * t87;
	t71 = 0.1e1 / t72 ^ 2;
	t99 = -t102 * t124 - t106 * t112;
	t117 = t102 * t130 + t108 * t99;
	t100 = -t102 * t125 + t106 * t115;
	t131 = t100 * t111;
	t89 = t117 * t114 - t131;
	t96 = t102 * t127 - t99 * t104;
	t83 = -t89 * t103 + t96 * t107;
	t136 = t71 * t83;
	t110 = sin(qJ(4));
	t119 = t103 * t96 + t107 * t89;
	t113 = cos(qJ(4));
	t90 = t100 * t114 + t117 * t111;
	t132 = t90 * t113;
	t77 = t119 * t110 + t132;
	t75 = 0.1e1 / t77 ^ 2;
	t133 = t90 * t110;
	t76 = -t119 * t113 + t133;
	t135 = t75 * t76;
	t86 = 0.1e1 / t87 ^ 2;
	t134 = t82 * t86;
	t129 = t104 * t107;
	t123 = t111 * t115;
	t122 = t112 * t114;
	t121 = t76 ^ 2 * t75 + 0.1e1;
	t120 = -t78 * t87 + t79 * t82;
	t91 = -t100 * t126 - t99 * t111;
	t118 = t100 * t103 * t104 + t107 * t91;
	t95 = -t111 * t128 + (-t108 * t123 - t122) * t105;
	t93 = (-(-t108 * t122 - t123) * t103 + t112 * t129) * t105;
	t92 = -t108 * t131 + t99 * t114;
	t88 = t116 * t111 - t98 * t114;
	t85 = 0.1e1 / t87;
	t84 = (-t97 * t111 - t98 * t126) * t103 - t98 * t129;
	t80 = 0.1e1 / (t82 ^ 2 * t86 + 0.1e1);
	t74 = 0.1e1 / t77;
	t73 = 0.1e1 / t121;
	t70 = 0.1e1 / t72;
	t69 = 0.1e1 / (t83 ^ 2 * t71 + 0.1e1);
	t68 = (t95 * t134 + t85 * t88) * t80 * t103;
	t67 = (-t93 * t134 + t84 * t85) * t80;
	t1 = [0, t67, t68, 0, 0, 0; 0, ((t100 * t129 - t91 * t103) * t70 - (t120 * t67 + t78 * t84 + t79 * t93) * t136) * t69, (t90 * t103 * t70 - (t120 * t68 + (t78 * t88 - t79 * t95) * t103) * t136) * t69, 0, 0, 0; 0, ((t92 * t110 - t118 * t113) * t74 - (t118 * t110 + t92 * t113) * t135) * t73, ((t107 * t132 + t89 * t110) * t74 - (-t107 * t133 + t89 * t113) * t135) * t73, t121 * t73, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:36
	% EndTime: 2019-10-09 23:23:37
	% DurationCPUTime: 0.81s
	% Computational Cost: add. (2578->75), mult. (7606->169), div. (95->9), fcn. (10284->19), ass. (0->89)
	t145 = sin(pkin(8));
	t146 = sin(pkin(7));
	t148 = cos(pkin(8));
	t149 = cos(pkin(7));
	t150 = cos(pkin(6));
	t154 = sin(qJ(2));
	t188 = cos(pkin(14));
	t168 = t188 * t154;
	t144 = sin(pkin(14));
	t158 = cos(qJ(2));
	t181 = t144 * t158;
	t140 = t150 * t168 + t181;
	t153 = sin(qJ(3));
	t157 = cos(qJ(3));
	t167 = t188 * t158;
	t182 = t144 * t154;
	t163 = -t150 * t167 + t182;
	t147 = sin(pkin(6));
	t169 = t147 * t188;
	t162 = -t146 * t169 - t163 * t149;
	t160 = -t140 * t153 + t162 * t157;
	t190 = t160 * t148 + (t163 * t146 - t149 * t169) * t145;
	t142 = -t150 * t182 + t167;
	t141 = -t150 * t181 - t168;
	t179 = t146 * t147;
	t164 = t141 * t149 + t144 * t179;
	t130 = -t142 * t153 + t164 * t157;
	t156 = cos(qJ(4));
	t177 = t148 * t156;
	t138 = t144 * t147 * t149 - t141 * t146;
	t183 = t138 * t145;
	t131 = t142 * t157 + t164 * t153;
	t152 = sin(qJ(4));
	t184 = t131 * t152;
	t113 = -t130 * t177 - t156 * t183 + t184;
	t129 = t140 * t157 + t162 * t153;
	t110 = t129 * t152 - t190 * t156;
	t172 = t154 * t157;
	t174 = t153 * t158;
	t178 = t146 * t150;
	t137 = t153 * t178 + (t149 * t174 + t172) * t147;
	t171 = t157 * t158;
	t173 = t154 * t153;
	t136 = t157 * t178 + (t149 * t171 - t173) * t147;
	t166 = t136 * t148 + (t150 * t149 - t158 * t179) * t145;
	t121 = t137 * t152 - t166 * t156;
	t109 = atan2(-t110, t121);
	t106 = sin(t109);
	t107 = cos(t109);
	t100 = -t106 * t110 + t107 * t121;
	t99 = 0.1e1 / t100 ^ 2;
	t189 = t113 * t99;
	t114 = t131 * t156 + (t130 * t148 + t183) * t152;
	t123 = -t130 * t145 + t138 * t148;
	t151 = sin(qJ(5));
	t155 = cos(qJ(5));
	t105 = t114 * t155 + t123 * t151;
	t103 = 0.1e1 / t105 ^ 2;
	t104 = t114 * t151 - t123 * t155;
	t187 = t103 * t104;
	t120 = 0.1e1 / t121 ^ 2;
	t186 = t110 * t120;
	t185 = t131 * t145;
	t180 = t146 * t145;
	t176 = t149 * t153;
	t175 = t149 * t157;
	t170 = t104 ^ 2 * t103 + 0.1e1;
	t132 = -t141 * t153 - t142 * t175;
	t165 = t132 * t148 + t142 * t180;
	t133 = t141 * t157 - t142 * t176;
	t126 = ((-t149 * t173 + t171) * t152 + (-(-t149 * t172 - t174) * t148 - t154 * t180) * t156) * t147;
	t125 = t142 * t146 * t148 - t132 * t145;
	t124 = t136 * t152 + t137 * t177;
	t122 = t137 * t156 + t166 * t152;
	t119 = 0.1e1 / t121;
	t118 = t130 * t156 - t148 * t184;
	t117 = t129 * t177 + t160 * t152;
	t116 = t133 * t156 + t165 * t152;
	t115 = (-t140 * t176 - t163 * t157) * t152 + (-(-t140 * t175 + t163 * t153) * t148 - t140 * t180) * t156;
	t112 = t129 * t156 + t190 * t152;
	t108 = 0.1e1 / (t110 ^ 2 * t120 + 0.1e1);
	t102 = 0.1e1 / t105;
	t101 = 0.1e1 / t170;
	t98 = 0.1e1 / t100;
	t97 = 0.1e1 / (t113 ^ 2 * t99 + 0.1e1);
	t96 = (-t115 * t119 + t126 * t186) * t108;
	t95 = (-t117 * t119 + t124 * t186) * t108;
	t94 = (-t112 * t119 + t122 * t186) * t108;
	t1 = [0, t96, t95, t94, 0, 0; 0, ((t133 * t152 - t165 * t156) * t98 - ((-t110 * t96 + t126) * t107 + (-t121 * t96 - t115) * t106) * t189) * t97, ((t130 * t152 + t131 * t177) * t98 - ((-t110 * t95 + t124) * t107 + (-t121 * t95 - t117) * t106) * t189) * t97, (t114 * t98 - ((-t110 * t94 + t122) * t107 + (-t121 * t94 - t112) * t106) * t189) * t97, 0, 0; 0, ((t116 * t151 - t125 * t155) * t102 - (t116 * t155 + t125 * t151) * t187) * t101, ((t118 * t151 - t155 * t185) * t102 - (t118 * t155 + t151 * t185) * t187) * t101, (-t151 * t102 + t155 * t187) * t113 * t101, t170 * t101, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:36
	% EndTime: 2019-10-09 23:23:37
	% DurationCPUTime: 1.31s
	% Computational Cost: add. (5587->99), mult. (16161->229), div. (125->9), fcn. (21882->21), ass. (0->109)
	t187 = sin(pkin(14));
	t191 = cos(pkin(14));
	t199 = sin(qJ(2));
	t194 = cos(pkin(6));
	t204 = cos(qJ(2));
	t222 = t194 * t204;
	t181 = -t187 * t199 + t191 * t222;
	t223 = t194 * t199;
	t182 = t187 * t204 + t191 * t223;
	t198 = sin(qJ(3));
	t203 = cos(qJ(3));
	t189 = sin(pkin(7));
	t190 = sin(pkin(6));
	t232 = t189 * t190;
	t217 = t191 * t232;
	t193 = cos(pkin(7));
	t225 = t193 * t198;
	t170 = t181 * t225 + t182 * t203 - t198 * t217;
	t197 = sin(qJ(4));
	t202 = cos(qJ(4));
	t188 = sin(pkin(8));
	t192 = cos(pkin(8));
	t207 = -t182 * t198 + (t181 * t193 - t217) * t203;
	t228 = t190 * t193;
	t212 = -t181 * t189 - t191 * t228;
	t205 = t212 * t188 + t207 * t192;
	t154 = t170 * t202 + t205 * t197;
	t196 = sin(qJ(5));
	t201 = cos(qJ(5));
	t206 = -t207 * t188 + t212 * t192;
	t140 = t154 * t196 - t206 * t201;
	t219 = t199 * t203;
	t220 = t198 * t204;
	t230 = t189 * t194;
	t178 = t198 * t230 + (t193 * t220 + t219) * t190;
	t218 = t203 * t204;
	t221 = t198 * t199;
	t177 = t203 * t230 + (t193 * t218 - t221) * t190;
	t180 = t194 * t193 - t204 * t232;
	t214 = t177 * t192 + t180 * t188;
	t166 = t178 * t202 + t214 * t197;
	t169 = -t177 * t188 + t180 * t192;
	t151 = t166 * t196 - t169 * t201;
	t139 = atan2(-t140, t151);
	t136 = sin(t139);
	t137 = cos(t139);
	t130 = -t136 * t140 + t137 * t151;
	t129 = 0.1e1 / t130 ^ 2;
	t184 = -t187 * t223 + t191 * t204;
	t183 = -t187 * t222 - t191 * t199;
	t210 = t183 * t193 + t187 * t232;
	t171 = -t184 * t198 + t210 * t203;
	t172 = t184 * t203 + t210 * t198;
	t211 = -t183 * t189 + t187 * t228;
	t209 = t211 * t188;
	t156 = t172 * t202 + (t171 * t192 + t209) * t197;
	t208 = -t171 * t188 + t211 * t192;
	t143 = t156 * t196 - t208 * t201;
	t238 = t129 * t143;
	t144 = t156 * t201 + t208 * t196;
	t226 = t192 * t202;
	t155 = -t171 * t226 + t172 * t197 - t202 * t209;
	t195 = sin(qJ(6));
	t200 = cos(qJ(6));
	t135 = t144 * t200 + t155 * t195;
	t133 = 0.1e1 / t135 ^ 2;
	t134 = t144 * t195 - t155 * t200;
	t237 = t133 * t134;
	t150 = 0.1e1 / t151 ^ 2;
	t236 = t140 * t150;
	t235 = t155 * t201;
	t234 = t188 * t201;
	t233 = t189 * t188;
	t231 = t189 * t192;
	t229 = t189 * t199;
	t227 = t192 * t197;
	t224 = t193 * t203;
	t216 = t134 ^ 2 * t133 + 0.1e1;
	t215 = -t136 * t151 - t137 * t140;
	t174 = -t183 * t198 - t184 * t224;
	t213 = t174 * t192 + t184 * t233;
	t175 = t183 * t203 - t184 * t225;
	t173 = -t181 * t198 - t182 * t224;
	t167 = -t174 * t188 + t184 * t231;
	t165 = -t178 * t197 + t214 * t202;
	t162 = (t177 * t202 - t178 * t227) * t196 - t178 * t234;
	t161 = ((t196 * t227 + t234) * (-t193 * t219 - t220) + ((-t193 * t221 + t218) * t202 + t188 * t197 * t229) * t196 - t192 * t201 * t229) * t190;
	t160 = t171 * t202 - t172 * t227;
	t159 = t171 * t197 + t172 * t226;
	t158 = t175 * t202 + t213 * t197;
	t157 = t175 * t197 - t213 * t202;
	t153 = -t170 * t197 + t205 * t202;
	t152 = t166 * t201 + t169 * t196;
	t149 = 0.1e1 / t151;
	t148 = t172 * t188 * t196 + t160 * t201;
	t147 = (-t170 * t227 + t207 * t202) * t196 - t170 * t234;
	t146 = t158 * t201 + t167 * t196;
	t145 = ((t181 * t203 - t182 * t225) * t202 + (t173 * t192 + t182 * t233) * t197) * t196 - (-t173 * t188 + t182 * t231) * t201;
	t142 = t154 * t201 + t206 * t196;
	t138 = 0.1e1 / (t140 ^ 2 * t150 + 0.1e1);
	t132 = 0.1e1 / t135;
	t131 = 0.1e1 / t216;
	t128 = 0.1e1 / t130;
	t127 = 0.1e1 / (t143 ^ 2 * t129 + 0.1e1);
	t126 = (-t149 * t153 + t165 * t236) * t196 * t138;
	t125 = (-t147 * t149 + t162 * t236) * t138;
	t124 = (-t145 * t149 + t161 * t236) * t138;
	t123 = (-t142 * t149 + t152 * t236) * t138;
	t1 = [0, t124, t125, t126, t123, 0; 0, ((t158 * t196 - t167 * t201) * t128 - (t215 * t124 - t136 * t145 + t137 * t161) * t238) * t127, ((t160 * t196 - t172 * t234) * t128 - (t215 * t125 - t136 * t147 + t137 * t162) * t238) * t127, (-t155 * t196 * t128 - ((-t136 * t153 + t137 * t165) * t196 + t215 * t126) * t238) * t127, (t144 * t128 - (t215 * t123 - t136 * t142 + t137 * t152) * t238) * t127, 0; 0, ((t146 * t195 - t157 * t200) * t132 - (t146 * t200 + t157 * t195) * t237) * t131, ((t148 * t195 - t159 * t200) * t132 - (t148 * t200 + t159 * t195) * t237) * t131, ((-t156 * t200 - t195 * t235) * t132 - (t156 * t195 - t200 * t235) * t237) * t131, (-t195 * t132 + t200 * t237) * t143 * t131, t216 * t131;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end