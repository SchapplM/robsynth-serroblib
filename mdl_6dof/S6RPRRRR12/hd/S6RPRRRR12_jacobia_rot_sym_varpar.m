% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR12
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
%   Wie in S6RPRRRR12_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:16
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR12_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR12_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobia_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (25->13), mult. (89->32), div. (20->9), fcn. (140->9), ass. (0->23)
	t28 = cos(pkin(6));
	t26 = sin(pkin(6));
	t30 = cos(qJ(1));
	t32 = t30 * t26;
	t22 = atan2(t32, t28);
	t19 = sin(t22);
	t20 = cos(t22);
	t14 = t19 * t32 + t20 * t28;
	t29 = sin(qJ(1));
	t37 = 0.1e1 / t14 ^ 2 * t29 ^ 2;
	t23 = t26 ^ 2;
	t21 = 0.1e1 / (0.1e1 + t30 ^ 2 * t23 / t28 ^ 2);
	t36 = t21 / t28;
	t25 = sin(pkin(14));
	t35 = t29 * t25;
	t27 = cos(pkin(14));
	t34 = t29 * t27;
	t33 = t30 * t25;
	t31 = t30 * t27;
	t18 = -t28 * t35 + t31;
	t17 = t28 * t34 + t33;
	t16 = 0.1e1 / t18 ^ 2;
	t1 = [-t29 * t26 * t36, 0, 0, 0, 0, 0; (0.1e1 / t14 * t32 - (-t20 * t23 * t30 * t36 + (t21 - 0.1e1) * t26 * t19) * t26 * t37) / (t23 * t37 + 0.1e1), 0, 0, 0, 0, 0; ((t28 * t31 - t35) / t18 - (-t28 * t33 - t34) * t17 * t16) / (t17 ^ 2 * t16 + 0.1e1), 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (156->21), mult. (450->48), div. (25->9), fcn. (637->13), ass. (0->38)
	t61 = cos(pkin(6));
	t59 = cos(pkin(14));
	t65 = cos(qJ(1));
	t69 = t65 * t59;
	t56 = sin(pkin(14));
	t63 = sin(qJ(1));
	t72 = t63 * t56;
	t51 = -t61 * t69 + t72;
	t57 = sin(pkin(7));
	t60 = cos(pkin(7));
	t58 = sin(pkin(6));
	t73 = t58 * t65;
	t46 = -t51 * t57 + t60 * t73;
	t50 = -t58 * t59 * t57 + t61 * t60;
	t45 = atan2(t46, t50);
	t42 = sin(t45);
	t43 = cos(t45);
	t37 = t42 * t46 + t43 * t50;
	t70 = t65 * t56;
	t71 = t63 * t59;
	t53 = -t61 * t71 - t70;
	t74 = t58 * t63;
	t47 = t53 * t57 - t60 * t74;
	t75 = t47 ^ 2 / t37 ^ 2;
	t54 = -t61 * t72 + t69;
	t62 = sin(qJ(3));
	t64 = cos(qJ(3));
	t66 = t53 * t60 + t57 * t74;
	t41 = t54 * t64 + t66 * t62;
	t39 = 0.1e1 / t41 ^ 2;
	t40 = t54 * t62 - t66 * t64;
	t68 = t40 ^ 2 * t39 + 0.1e1;
	t67 = t51 * t60 + t57 * t73;
	t52 = -t61 * t70 - t71;
	t49 = 0.1e1 / t50;
	t44 = 0.1e1 / (0.1e1 + t46 ^ 2 / t50 ^ 2);
	t38 = 0.1e1 / t68;
	t1 = [t47 * t49 * t44, 0, 0, 0, 0, 0; (t46 / t37 + (t42 + (t43 * t46 * t49 - t42) * t44) * t75) / (0.1e1 + t75), 0, 0, 0, 0, 0; ((t52 * t62 - t67 * t64) / t41 - (t52 * t64 + t67 * t62) * t40 * t39) * t38, 0, t68 * t38, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:47
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (950->43), mult. (2744->104), div. (55->9), fcn. (3759->17), ass. (0->63)
	t102 = cos(pkin(8));
	t103 = cos(pkin(7));
	t106 = sin(qJ(3));
	t109 = cos(qJ(3));
	t100 = sin(pkin(6));
	t110 = cos(qJ(1));
	t122 = t100 * t110;
	t99 = sin(pkin(7));
	t118 = t99 * t122;
	t104 = cos(pkin(6));
	t101 = cos(pkin(14));
	t119 = t110 * t101;
	t107 = sin(qJ(1));
	t97 = sin(pkin(14));
	t129 = t107 * t97;
	t92 = -t104 * t119 + t129;
	t120 = t107 * t101;
	t128 = t110 * t97;
	t93 = t104 * t128 + t120;
	t84 = (t103 * t92 + t118) * t109 + t93 * t106;
	t91 = t103 * t122 - t92 * t99;
	t98 = sin(pkin(8));
	t77 = t91 * t102 - t84 * t98;
	t105 = sin(qJ(4));
	t123 = t100 * t107;
	t94 = -t104 * t120 - t128;
	t114 = t103 * t123 - t94 * t99;
	t115 = t103 * t94 + t99 * t123;
	t95 = -t104 * t129 + t119;
	t86 = -t95 * t106 + t115 * t109;
	t111 = t86 * t102 + t114 * t98;
	t108 = cos(qJ(4));
	t87 = t115 * t106 + t95 * t109;
	t126 = t87 * t108;
	t72 = t111 * t105 + t126;
	t70 = 0.1e1 / t72 ^ 2;
	t127 = t87 * t105;
	t71 = -t111 * t108 + t127;
	t133 = t70 * t71;
	t130 = t104 * t99;
	t82 = -(t109 * t130 + (t101 * t103 * t109 - t106 * t97) * t100) * t98 + (-t100 * t101 * t99 + t104 * t103) * t102;
	t76 = atan2(t77, t82);
	t74 = cos(t76);
	t132 = t74 * t77;
	t73 = sin(t76);
	t67 = t73 * t77 + t74 * t82;
	t66 = 0.1e1 / t67 ^ 2;
	t78 = -t114 * t102 + t86 * t98;
	t131 = t78 ^ 2 * t66;
	t121 = t103 * t106;
	t117 = t71 ^ 2 * t70 + 0.1e1;
	t116 = t102 * t84 + t91 * t98;
	t112 = t106 * t118 - t93 * t109 + t92 * t121;
	t90 = -t106 * t130 + (-t101 * t121 - t109 * t97) * t100;
	t81 = 0.1e1 / t82 ^ 2;
	t80 = 0.1e1 / t82;
	t75 = 0.1e1 / (t77 ^ 2 * t81 + 0.1e1);
	t69 = 0.1e1 / t72;
	t68 = 0.1e1 / t117;
	t65 = 0.1e1 / t67;
	t64 = 0.1e1 / (0.1e1 + t131);
	t63 = (t77 * t81 * t90 + t112 * t80) * t98 * t75;
	t1 = [t78 * t80 * t75, 0, t63, 0, 0, 0; (t77 * t65 + (t73 + (t80 * t132 - t73) * t75) * t131) * t64, 0, (t87 * t98 * t65 + ((t112 * t73 - t74 * t90) * t98 + (-t73 * t82 + t132) * t63) * t78 * t66) * t64, 0, 0, 0; ((t105 * t112 - t116 * t108) * t69 - (t116 * t105 + t108 * t112) * t133) * t68, 0, ((t102 * t126 + t86 * t105) * t69 - (-t102 * t127 + t86 * t108) * t133) * t68, t117 * t68, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:47
	% DurationCPUTime: 0.96s
	% Computational Cost: add. (2433->62), mult. (7147->139), div. (85->9), fcn. (9679->19), ass. (0->78)
	t140 = cos(pkin(14));
	t136 = sin(pkin(14));
	t147 = sin(qJ(1));
	t168 = t147 * t136;
	t143 = cos(pkin(6));
	t151 = cos(qJ(1));
	t169 = t143 * t151;
	t132 = -t140 * t169 + t168;
	t166 = t147 * t140;
	t133 = t136 * t169 + t166;
	t146 = sin(qJ(3));
	t150 = cos(qJ(3));
	t138 = sin(pkin(7));
	t139 = sin(pkin(6));
	t172 = t139 * t151;
	t164 = t138 * t172;
	t142 = cos(pkin(7));
	t170 = t142 * t146;
	t124 = t132 * t170 - t133 * t150 + t146 * t164;
	t145 = sin(qJ(4));
	t149 = cos(qJ(4));
	t123 = (t132 * t142 + t164) * t150 + t133 * t146;
	t130 = -t132 * t138 + t142 * t172;
	t137 = sin(pkin(8));
	t141 = cos(pkin(8));
	t163 = t123 * t141 + t130 * t137;
	t108 = t124 * t149 + t163 * t145;
	t105 = -t124 * t145 + t163 * t149;
	t134 = t140 * t151 - t143 * t168;
	t160 = -t136 * t151 - t143 * t166;
	t167 = t147 * t139;
	t156 = t138 * t167 + t160 * t142;
	t125 = t134 * t150 + t156 * t146;
	t153 = t134 * t146 - t156 * t150;
	t157 = -t160 * t138 + t142 * t167;
	t154 = t157 * t137;
	t110 = t125 * t149 + (-t153 * t141 + t154) * t145;
	t118 = t153 * t137 + t157 * t141;
	t144 = sin(qJ(5));
	t148 = cos(qJ(5));
	t100 = t110 * t148 + t118 * t144;
	t98 = 0.1e1 / t100 ^ 2;
	t99 = t110 * t144 - t118 * t148;
	t181 = t98 * t99;
	t152 = t153 * t149;
	t176 = t125 * t145;
	t109 = t141 * t152 - t149 * t154 + t176;
	t173 = t138 * t143;
	t129 = t146 * t173 + (t136 * t150 + t140 * t170) * t139;
	t128 = t150 * t173 + (t140 * t142 * t150 - t136 * t146) * t139;
	t162 = t128 * t141 + (-t138 * t139 * t140 + t142 * t143) * t137;
	t115 = t129 * t145 - t162 * t149;
	t104 = atan2(-t105, t115);
	t101 = sin(t104);
	t102 = cos(t104);
	t95 = -t101 * t105 + t102 * t115;
	t94 = 0.1e1 / t95 ^ 2;
	t180 = t109 * t94;
	t179 = t109 ^ 2 * t94;
	t114 = 0.1e1 / t115 ^ 2;
	t178 = t105 * t114;
	t177 = t125 * t137;
	t171 = t141 * t149;
	t165 = t98 * t99 ^ 2 + 0.1e1;
	t119 = t128 * t145 + t129 * t171;
	t117 = -t123 * t137 + t130 * t141;
	t116 = t129 * t149 + t162 * t145;
	t113 = 0.1e1 / t115;
	t112 = -t141 * t176 - t152;
	t111 = -t123 * t145 - t124 * t171;
	t103 = 0.1e1 / (t105 ^ 2 * t114 + 0.1e1);
	t97 = 0.1e1 / t100;
	t96 = 0.1e1 / t165;
	t93 = 0.1e1 / t95;
	t92 = 0.1e1 / (0.1e1 + t179);
	t91 = (-t111 * t113 + t119 * t178) * t103;
	t90 = (t108 * t113 + t116 * t178) * t103;
	t1 = [-t109 * t113 * t103, 0, t91, t90, 0, 0; (-t105 * t93 - (-t101 + (t102 * t105 * t113 + t101) * t103) * t179) * t92, 0, ((t125 * t171 - t153 * t145) * t93 - ((-t105 * t91 + t119) * t102 + (-t115 * t91 - t111) * t101) * t180) * t92, (t110 * t93 - ((-t105 * t90 + t116) * t102 + (-t115 * t90 + t108) * t101) * t180) * t92, 0, 0; ((t108 * t144 - t117 * t148) * t97 - (t108 * t148 + t117 * t144) * t181) * t96, 0, ((t112 * t144 - t148 * t177) * t97 - (t112 * t148 + t144 * t177) * t181) * t96, (-t144 * t97 + t148 * t181) * t96 * t109, t165 * t96, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:47
	% EndTime: 2019-10-10 09:15:49
	% DurationCPUTime: 1.78s
	% Computational Cost: add. (5345->78), mult. (15419->180), div. (115->9), fcn. (20897->21), ass. (0->95)
	t187 = cos(pkin(6));
	t180 = sin(pkin(14));
	t234 = sin(qJ(1));
	t214 = t234 * t180;
	t184 = cos(pkin(14));
	t196 = cos(qJ(1));
	t217 = t196 * t184;
	t176 = -t187 * t217 + t214;
	t213 = t234 * t184;
	t218 = t196 * t180;
	t177 = t187 * t218 + t213;
	t191 = sin(qJ(3));
	t195 = cos(qJ(3));
	t182 = sin(pkin(7));
	t183 = sin(pkin(6));
	t221 = t183 * t196;
	t216 = t182 * t221;
	t186 = cos(pkin(7));
	t219 = t186 * t191;
	t169 = t176 * t219 - t177 * t195 + t191 * t216;
	t190 = sin(qJ(4));
	t194 = cos(qJ(4));
	t168 = (t176 * t186 + t216) * t195 + t177 * t191;
	t174 = -t176 * t182 + t186 * t221;
	t181 = sin(pkin(8));
	t185 = cos(pkin(8));
	t210 = t168 * t185 + t174 * t181;
	t154 = t169 * t194 + t210 * t190;
	t189 = sin(qJ(5));
	t193 = cos(qJ(5));
	t201 = t168 * t181 - t174 * t185;
	t240 = t154 * t189 + t201 * t193;
	t142 = t154 * t193 - t201 * t189;
	t237 = t169 * t190 - t210 * t194;
	t222 = t182 * t187;
	t173 = t191 * t222 + (t180 * t195 + t184 * t219) * t183;
	t172 = t195 * t222 + (t184 * t186 * t195 - t180 * t191) * t183;
	t175 = -t183 * t184 * t182 + t187 * t186;
	t209 = t172 * t185 + t175 * t181;
	t161 = t173 * t194 + t209 * t190;
	t166 = -t172 * t181 + t175 * t185;
	t149 = t161 * t189 - t166 * t193;
	t138 = atan2(t240, t149);
	t133 = sin(t138);
	t134 = cos(t138);
	t129 = t133 * t240 + t134 * t149;
	t128 = 0.1e1 / t129 ^ 2;
	t207 = -t187 * t213 - t218;
	t215 = t183 * t234;
	t203 = t182 * t215 + t207 * t186;
	t206 = -t187 * t214 + t217;
	t199 = t206 * t191 - t203 * t195;
	t204 = -t207 * t182 + t186 * t215;
	t202 = t204 * t181;
	t170 = t203 * t191 + t206 * t195;
	t227 = t170 * t194;
	t156 = t227 + (-t199 * t185 + t202) * t190;
	t197 = t199 * t181 + t204 * t185;
	t143 = t156 * t189 - t197 * t193;
	t233 = t128 * t143;
	t144 = t156 * t193 + t197 * t189;
	t198 = t199 * t194;
	t155 = t170 * t190 + t185 * t198 - t194 * t202;
	t188 = sin(qJ(6));
	t192 = cos(qJ(6));
	t136 = t144 * t192 + t155 * t188;
	t132 = 0.1e1 / t136 ^ 2;
	t135 = t144 * t188 - t155 * t192;
	t232 = t132 * t135;
	t231 = t134 * t240;
	t148 = 0.1e1 / t149 ^ 2;
	t230 = t240 * t148;
	t229 = t143 ^ 2 * t128;
	t228 = t155 * t193;
	t223 = t181 * t193;
	t220 = t185 * t190;
	t212 = t135 ^ 2 * t132 + 0.1e1;
	t211 = -t133 * t149 + t231;
	t160 = -t173 * t190 + t209 * t194;
	t159 = (t172 * t194 - t173 * t220) * t189 - t173 * t223;
	t158 = -t170 * t220 - t198;
	t157 = t185 * t227 - t199 * t190;
	t150 = t161 * t193 + t166 * t189;
	t147 = 0.1e1 / t149;
	t146 = t170 * t181 * t189 + t158 * t193;
	t145 = (-t168 * t194 + t169 * t220) * t189 + t169 * t223;
	t137 = 0.1e1 / (t148 * t240 ^ 2 + 0.1e1);
	t131 = 0.1e1 / t136;
	t130 = 0.1e1 / t212;
	t127 = 0.1e1 / t129;
	t126 = 0.1e1 / (0.1e1 + t229);
	t125 = (-t147 * t237 - t160 * t230) * t189 * t137;
	t124 = (-t145 * t147 - t159 * t230) * t137;
	t123 = (t142 * t147 - t150 * t230) * t137;
	t1 = [-t143 * t147 * t137, 0, t124, t125, t123, 0; (t240 * t127 - (-t133 + (-t147 * t231 + t133) * t137) * t229) * t126, 0, ((t158 * t189 - t170 * t223) * t127 - (t211 * t124 - t133 * t145 + t134 * t159) * t233) * t126, (-t155 * t189 * t127 - ((-t133 * t237 + t134 * t160) * t189 + t211 * t125) * t233) * t126, (t144 * t127 - (t211 * t123 + t133 * t142 + t134 * t150) * t233) * t126, 0; ((t142 * t188 - t192 * t237) * t131 - (t142 * t192 + t188 * t237) * t232) * t130, 0, ((t146 * t188 - t157 * t192) * t131 - (t146 * t192 + t157 * t188) * t232) * t130, ((-t156 * t192 - t188 * t228) * t131 - (t156 * t188 - t192 * t228) * t232) * t130, (-t188 * t131 + t192 * t232) * t143 * t130, t212 * t130;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end