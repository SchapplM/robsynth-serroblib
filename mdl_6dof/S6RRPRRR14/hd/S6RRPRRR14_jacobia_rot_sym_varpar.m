% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR14
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
%   Wie in S6RRPRRR14_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR14_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR14_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (35->13), mult. (116->34), div. (25->9), fcn. (175->9), ass. (0->25)
	t30 = cos(pkin(6));
	t29 = sin(pkin(6));
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
	t1 = [-t32 * t29 * t41, 0, 0, 0, 0, 0; (0.1e1 / t18 * t38 - (-t24 * t27 * t34 * t41 + (t25 - 0.1e1) * t29 * t23) * t29 * t42) / (t27 * t42 + 0.1e1), 0, 0, 0, 0, 0; ((t30 * t36 - t40) / t22 - (-t30 * t37 - t39) * t21 * t20) * t19, t35 * t19, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (273->29), mult. (795->75), div. (50->9), fcn. (1135->13), ass. (0->47)
	t61 = cos(pkin(6));
	t64 = cos(qJ(2));
	t65 = cos(qJ(1));
	t68 = t64 * t65;
	t62 = sin(qJ(2));
	t63 = sin(qJ(1));
	t70 = t63 * t62;
	t50 = -t61 * t68 + t70;
	t57 = sin(pkin(7));
	t60 = cos(pkin(7));
	t58 = sin(pkin(6));
	t72 = t58 * t65;
	t44 = -t50 * t57 + t60 * t72;
	t49 = -t57 * t58 * t64 + t60 * t61;
	t43 = atan2(t44, t49);
	t40 = sin(t43);
	t41 = cos(t43);
	t34 = t40 * t44 + t41 * t49;
	t33 = 0.1e1 / t34 ^ 2;
	t69 = t63 * t64;
	t71 = t62 * t65;
	t52 = -t61 * t69 - t71;
	t73 = t58 * t63;
	t45 = t52 * t57 - t60 * t73;
	t78 = t33 * t45 ^ 2;
	t53 = -t61 * t70 + t68;
	t56 = sin(pkin(14));
	t59 = cos(pkin(14));
	t66 = t52 * t60 + t57 * t73;
	t39 = t53 * t59 + t66 * t56;
	t37 = 0.1e1 / t39 ^ 2;
	t38 = t53 * t56 - t66 * t59;
	t77 = t37 * t38;
	t76 = t41 * t44;
	t75 = t53 * t60;
	t74 = t58 * t62;
	t67 = t50 * t60 + t57 * t72;
	t51 = -t61 * t71 - t69;
	t48 = 0.1e1 / t49 ^ 2;
	t47 = 0.1e1 / t49;
	t42 = 0.1e1 / (t44 ^ 2 * t48 + 0.1e1);
	t36 = 0.1e1 / t39;
	t35 = 0.1e1 / (t37 * t38 ^ 2 + 0.1e1);
	t32 = 0.1e1 / t34;
	t31 = 0.1e1 / (0.1e1 + t78);
	t30 = (-t44 * t48 * t74 + t47 * t51) * t57 * t42;
	t1 = [t45 * t47 * t42, t30, 0, 0, 0, 0; (t44 * t32 + (t40 + (t47 * t76 - t40) * t42) * t78) * t31, (t53 * t57 * t32 + ((t40 * t51 + t41 * t74) * t57 + (-t40 * t49 + t76) * t30) * t45 * t33) * t31, 0, 0, 0, 0; ((t51 * t56 - t67 * t59) * t36 - (t51 * t59 + t67 * t56) * t77) * t35, ((t52 * t56 + t59 * t75) * t36 - (t52 * t59 - t56 * t75) * t77) * t35, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:59
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (955->48), mult. (2758->114), div. (55->9), fcn. (3778->17), ass. (0->66)
	t117 = cos(pkin(6));
	t122 = cos(qJ(2));
	t123 = cos(qJ(1));
	t131 = t123 * t122;
	t119 = sin(qJ(2));
	t120 = sin(qJ(1));
	t134 = t120 * t119;
	t105 = -t117 * t131 + t134;
	t112 = sin(pkin(7));
	t116 = cos(pkin(7));
	t113 = sin(pkin(6));
	t136 = t113 * t123;
	t104 = -t105 * t112 + t116 * t136;
	t111 = sin(pkin(8));
	t115 = cos(pkin(8));
	t132 = t123 * t119;
	t133 = t120 * t122;
	t106 = t117 * t132 + t133;
	t110 = sin(pkin(14));
	t114 = cos(pkin(14));
	t127 = t105 * t116 + t112 * t136;
	t95 = t106 * t110 + t127 * t114;
	t88 = t104 * t115 - t95 * t111;
	t118 = sin(qJ(4));
	t121 = cos(qJ(4));
	t107 = -t117 * t133 - t132;
	t137 = t113 * t120;
	t126 = -t107 * t112 + t116 * t137;
	t125 = t107 * t116 + t112 * t137;
	t108 = -t117 * t134 + t131;
	t139 = t108 * t110;
	t97 = t125 * t114 - t139;
	t124 = t126 * t111 + t97 * t115;
	t98 = t108 * t114 + t125 * t110;
	t83 = t124 * t118 + t98 * t121;
	t81 = 0.1e1 / t83 ^ 2;
	t82 = t98 * t118 - t124 * t121;
	t144 = t81 * t82;
	t135 = t114 * t116;
	t94 = -(t117 * t112 * t114 + (-t110 * t119 + t122 * t135) * t113) * t111 + (-t113 * t122 * t112 + t117 * t116) * t115;
	t87 = atan2(t88, t94);
	t85 = cos(t87);
	t143 = t85 * t88;
	t84 = sin(t87);
	t78 = t84 * t88 + t85 * t94;
	t77 = 0.1e1 / t78 ^ 2;
	t89 = t97 * t111 - t126 * t115;
	t142 = t89 ^ 2 * t77;
	t138 = t112 * t115;
	t130 = t81 * t82 ^ 2 + 0.1e1;
	t129 = t104 * t111 + t115 * t95;
	t99 = -t107 * t110 - t108 * t135;
	t128 = t108 * t111 * t112 + t115 * t99;
	t101 = (-(-t110 * t122 - t119 * t135) * t111 + t119 * t138) * t113;
	t100 = t107 * t114 - t116 * t139;
	t96 = -t106 * t114 + t127 * t110;
	t93 = 0.1e1 / t94 ^ 2;
	t92 = 0.1e1 / t94;
	t91 = (t105 * t110 - t106 * t135) * t111 - t106 * t138;
	t86 = 0.1e1 / (t88 ^ 2 * t93 + 0.1e1);
	t80 = 0.1e1 / t83;
	t79 = 0.1e1 / t130;
	t76 = 0.1e1 / t78;
	t75 = 0.1e1 / (0.1e1 + t142);
	t74 = (-t101 * t88 * t93 + t91 * t92) * t86;
	t1 = [t89 * t92 * t86, t74, 0, 0, 0, 0; (t88 * t76 + (t84 + (t92 * t143 - t84) * t86) * t142) * t75, ((t108 * t138 - t99 * t111) * t76 + (t85 * t101 + t84 * t91 + (-t84 * t94 + t143) * t74) * t89 * t77) * t75, 0, 0, 0, 0; ((t96 * t118 - t129 * t121) * t80 - (t129 * t118 + t96 * t121) * t144) * t79, ((t100 * t118 - t128 * t121) * t80 - (t100 * t121 + t128 * t118) * t144) * t79, 0, t130 * t79, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:59
	% EndTime: 2019-10-10 11:11:00
	% DurationCPUTime: 1.04s
	% Computational Cost: add. (2431->68), mult. (7151->158), div. (85->9), fcn. (9681->19), ass. (0->84)
	t161 = cos(pkin(6));
	t167 = cos(qJ(2));
	t203 = sin(qJ(1));
	t183 = t203 * t167;
	t164 = sin(qJ(2));
	t168 = cos(qJ(1));
	t187 = t168 * t164;
	t151 = t161 * t187 + t183;
	t154 = sin(pkin(14));
	t158 = cos(pkin(14));
	t184 = t203 * t164;
	t186 = t168 * t167;
	t150 = -t161 * t186 + t184;
	t156 = sin(pkin(7));
	t160 = cos(pkin(7));
	t157 = sin(pkin(6));
	t191 = t157 * t168;
	t177 = t150 * t160 + t156 * t191;
	t140 = -t151 * t158 + t177 * t154;
	t163 = sin(qJ(4));
	t166 = cos(qJ(4));
	t139 = t151 * t154 + t177 * t158;
	t148 = -t150 * t156 + t160 * t191;
	t155 = sin(pkin(8));
	t159 = cos(pkin(8));
	t204 = -t139 * t159 - t148 * t155;
	t122 = -t140 * t166 + t204 * t163;
	t209 = t140 * t163 + t204 * t166;
	t152 = -t161 * t184 + t186;
	t176 = t161 * t183 + t187;
	t185 = t157 * t203;
	t172 = t156 * t185 - t176 * t160;
	t170 = t152 * t154 - t172 * t158;
	t173 = t176 * t156 + t160 * t185;
	t205 = -t173 * t155 + t170 * t159;
	t189 = t160 * t167;
	t192 = t156 * t161;
	t147 = t157 * t164 * t158 + (t157 * t189 + t192) * t154;
	t179 = (t158 * t192 + (-t154 * t164 + t158 * t189) * t157) * t159 + (-t157 * t167 * t156 + t161 * t160) * t155;
	t130 = t147 * t163 - t179 * t166;
	t119 = atan2(t209, t130);
	t112 = sin(t119);
	t113 = cos(t119);
	t110 = t112 * t209 + t113 * t130;
	t109 = 0.1e1 / t110 ^ 2;
	t141 = t152 * t158 + t172 * t154;
	t124 = t141 * t163 + t205 * t166;
	t202 = t109 * t124;
	t201 = t113 * t209;
	t125 = t141 * t166 - t205 * t163;
	t133 = t170 * t155 + t173 * t159;
	t162 = sin(qJ(5));
	t165 = cos(qJ(5));
	t118 = t125 * t165 + t133 * t162;
	t115 = 0.1e1 / t118 ^ 2;
	t117 = t125 * t162 - t133 * t165;
	t200 = t115 * t117;
	t129 = 0.1e1 / t130 ^ 2;
	t199 = t209 * t129;
	t198 = t124 ^ 2 * t109;
	t194 = t154 * t160;
	t193 = t156 * t155;
	t190 = t158 * t160;
	t188 = t164 * t160;
	t182 = t117 ^ 2 * t115 + 0.1e1;
	t181 = -t112 * t130 + t201;
	t142 = -t152 * t190 + t176 * t154;
	t178 = t142 * t159 + t152 * t193;
	t143 = -t152 * t194 - t176 * t158;
	t135 = t152 * t156 * t159 - t142 * t155;
	t134 = ((-t154 * t188 + t158 * t167) * t163 + (-(-t154 * t167 - t158 * t188) * t159 - t164 * t193) * t166) * t157;
	t132 = -t139 * t155 + t148 * t159;
	t131 = t147 * t166 + t179 * t163;
	t128 = 0.1e1 / t130;
	t127 = t143 * t166 + t178 * t163;
	t126 = (-t150 * t158 - t151 * t194) * t163 + (-(t150 * t154 - t151 * t190) * t159 - t151 * t193) * t166;
	t116 = 0.1e1 / (t129 * t209 ^ 2 + 0.1e1);
	t114 = 0.1e1 / t118;
	t111 = 0.1e1 / t182;
	t108 = 0.1e1 / t110;
	t107 = 0.1e1 / (0.1e1 + t198);
	t106 = (-t126 * t128 - t134 * t199) * t116;
	t105 = (-t122 * t128 - t131 * t199) * t116;
	t1 = [-t124 * t128 * t116, t106, 0, t105, 0, 0; (t209 * t108 - (-t112 + (-t128 * t201 + t112) * t116) * t198) * t107, ((t143 * t163 - t178 * t166) * t108 - (t181 * t106 - t112 * t126 + t113 * t134) * t202) * t107, 0, (t125 * t108 - (t181 * t105 - t112 * t122 + t113 * t131) * t202) * t107, 0, 0; ((-t122 * t162 - t132 * t165) * t114 - (-t122 * t165 + t132 * t162) * t200) * t111, ((t127 * t162 - t135 * t165) * t114 - (t127 * t165 + t135 * t162) * t200) * t111, 0, (-t162 * t114 + t165 * t200) * t124 * t111, t182 * t111, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:59
	% EndTime: 2019-10-10 11:11:01
	% DurationCPUTime: 1.98s
	% Computational Cost: add. (5348->88), mult. (15437->209), div. (115->9), fcn. (20918->21), ass. (0->100)
	t203 = sin(qJ(2));
	t207 = cos(qJ(2));
	t208 = cos(qJ(1));
	t247 = cos(pkin(6));
	t226 = t208 * t247;
	t248 = sin(qJ(1));
	t189 = t203 * t226 + t248 * t207;
	t193 = sin(pkin(14));
	t197 = cos(pkin(14));
	t188 = t248 * t203 - t207 * t226;
	t195 = sin(pkin(7));
	t199 = cos(pkin(7));
	t196 = sin(pkin(6));
	t232 = t196 * t208;
	t219 = t188 * t199 + t195 * t232;
	t177 = -t189 * t197 + t219 * t193;
	t202 = sin(qJ(4));
	t206 = cos(qJ(4));
	t176 = t189 * t193 + t219 * t197;
	t185 = -t188 * t195 + t199 * t232;
	t194 = sin(pkin(8));
	t198 = cos(pkin(8));
	t222 = t176 * t198 + t185 * t194;
	t161 = t177 * t206 + t222 * t202;
	t201 = sin(qJ(5));
	t205 = cos(qJ(5));
	t212 = t176 * t194 - t185 * t198;
	t255 = t161 * t201 + t212 * t205;
	t149 = t161 * t205 - t212 * t201;
	t252 = t177 * t202 - t222 * t206;
	t224 = t247 * t248;
	t190 = -t203 * t224 + t208 * t207;
	t218 = t208 * t203 + t207 * t224;
	t228 = t196 * t248;
	t215 = t195 * t228 - t218 * t199;
	t213 = t190 * t193 - t215 * t197;
	t216 = t218 * t195 + t199 * t228;
	t249 = -t216 * t194 + t213 * t198;
	t225 = t247 * t195;
	t229 = t199 * t207;
	t184 = t196 * t203 * t197 + (t196 * t229 + t225) * t193;
	t183 = t197 * t225 + (-t193 * t203 + t197 * t229) * t196;
	t187 = -t196 * t207 * t195 + t247 * t199;
	t221 = t183 * t198 + t187 * t194;
	t168 = t184 * t206 + t221 * t202;
	t174 = -t183 * t194 + t187 * t198;
	t156 = t168 * t201 - t174 * t205;
	t143 = atan2(t255, t156);
	t138 = sin(t143);
	t139 = cos(t143);
	t136 = t138 * t255 + t139 * t156;
	t135 = 0.1e1 / t136 ^ 2;
	t178 = t190 * t197 + t215 * t193;
	t163 = t178 * t206 - t249 * t202;
	t209 = t213 * t194 + t216 * t198;
	t150 = t163 * t201 - t209 * t205;
	t246 = t135 * t150;
	t245 = t139 * t255;
	t151 = t163 * t205 + t209 * t201;
	t162 = t178 * t202 + t249 * t206;
	t200 = sin(qJ(6));
	t204 = cos(qJ(6));
	t145 = t151 * t204 + t162 * t200;
	t142 = 0.1e1 / t145 ^ 2;
	t144 = t151 * t200 - t162 * t204;
	t244 = t142 * t144;
	t155 = 0.1e1 / t156 ^ 2;
	t243 = t255 * t155;
	t242 = t150 ^ 2 * t135;
	t241 = t162 * t205;
	t236 = t193 * t199;
	t235 = t195 * t194;
	t234 = t195 * t198;
	t233 = t195 * t203;
	t231 = t197 * t199;
	t230 = t199 * t203;
	t227 = t144 ^ 2 * t142 + 0.1e1;
	t223 = -t138 * t156 + t245;
	t180 = -t190 * t231 + t218 * t193;
	t220 = t180 * t198 + t190 * t235;
	t181 = -t190 * t236 - t218 * t197;
	t179 = t188 * t193 - t189 * t231;
	t172 = -t180 * t194 + t190 * t234;
	t167 = -t184 * t202 + t221 * t206;
	t166 = ((t198 * t202 * t201 + t194 * t205) * (-t193 * t207 - t197 * t230) + ((-t193 * t230 + t197 * t207) * t206 + t194 * t202 * t233) * t201 - t198 * t205 * t233) * t196;
	t165 = t181 * t206 + t220 * t202;
	t164 = t181 * t202 - t220 * t206;
	t157 = t168 * t205 + t174 * t201;
	t154 = 0.1e1 / t156;
	t153 = t165 * t205 + t172 * t201;
	t152 = ((-t188 * t197 - t189 * t236) * t206 + (t179 * t198 + t189 * t235) * t202) * t201 - (-t179 * t194 + t189 * t234) * t205;
	t141 = 0.1e1 / t145;
	t140 = 0.1e1 / (t155 * t255 ^ 2 + 0.1e1);
	t137 = 0.1e1 / t227;
	t134 = 0.1e1 / t136;
	t133 = 0.1e1 / (0.1e1 + t242);
	t132 = (-t154 * t252 - t167 * t243) * t201 * t140;
	t131 = (-t152 * t154 - t166 * t243) * t140;
	t130 = (t149 * t154 - t157 * t243) * t140;
	t1 = [-t150 * t154 * t140, t131, 0, t132, t130, 0; (t255 * t134 - (-t138 + (-t154 * t245 + t138) * t140) * t242) * t133, ((t165 * t201 - t172 * t205) * t134 - (t223 * t131 - t138 * t152 + t139 * t166) * t246) * t133, 0, (-t162 * t201 * t134 - ((-t138 * t252 + t139 * t167) * t201 + t223 * t132) * t246) * t133, (t151 * t134 - (t223 * t130 + t138 * t149 + t139 * t157) * t246) * t133, 0; ((t149 * t200 - t204 * t252) * t141 - (t149 * t204 + t200 * t252) * t244) * t137, ((t153 * t200 - t164 * t204) * t141 - (t153 * t204 + t164 * t200) * t244) * t137, 0, ((-t163 * t204 - t200 * t241) * t141 - (t163 * t200 - t204 * t241) * t244) * t137, (-t200 * t141 + t204 * t244) * t150 * t137, t227 * t137;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end