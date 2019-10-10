% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR7
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
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RRRRPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:42
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:15
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (215->39), mult. (853->106), div. (126->12), fcn. (1047->9), ass. (0->54)
	t99 = sin(pkin(6));
	t93 = t99 ^ 2;
	t100 = cos(pkin(6));
	t95 = 0.1e1 / t100 ^ 2;
	t104 = cos(qJ(1));
	t98 = t104 ^ 2;
	t89 = t98 * t93 * t95 + 0.1e1;
	t102 = sin(qJ(1));
	t97 = t102 ^ 2;
	t126 = 0.1e1 / t89 ^ 2 * t97;
	t131 = t126 * t95;
	t122 = t104 * t99;
	t88 = atan2(t122, t100);
	t84 = sin(t88);
	t85 = cos(t88);
	t72 = t85 * t100 + t84 * t122;
	t67 = 0.1e1 / t72;
	t103 = cos(qJ(2));
	t118 = t104 * t103;
	t101 = sin(qJ(2));
	t121 = t102 * t101;
	t113 = t100 * t121 - t118;
	t77 = 0.1e1 / t113;
	t94 = 0.1e1 / t100;
	t68 = 0.1e1 / t72 ^ 2;
	t78 = 0.1e1 / t113 ^ 2;
	t119 = t104 * t101;
	t120 = t102 * t103;
	t81 = -t100 * t119 - t120;
	t82 = t100 * t120 + t119;
	t71 = t81 * qJD(1) - t82 * qJD(2);
	t128 = t71 * t77 * t78;
	t115 = t100 * t118;
	t70 = -qJD(1) * t115 - qJD(2) * t118 + (qJD(2) * t100 + qJD(1)) * t121;
	t129 = t70 * t78;
	t76 = t82 ^ 2;
	t75 = t76 * t78 + 0.1e1;
	t130 = (t76 * t128 - t82 * t129) / t75 ^ 2;
	t127 = t81 * t82;
	t125 = t93 * t94;
	t124 = t102 * t68;
	t123 = t104 * t68;
	t117 = qJD(1) * t104;
	t86 = 0.1e1 / t89;
	t116 = (t86 - 0.1e1) * t99;
	t114 = -0.2e1 * t94 * t131;
	t80 = t115 - t121;
	t63 = (-t104 * t85 * t86 * t125 + t84 * t116) * t102;
	t92 = t99 * t93;
	t73 = 0.1e1 / t75;
	t69 = t67 * t68;
	t66 = t97 * t93 * t68 + 0.1e1;
	t62 = qJD(1) * t63;
	t1 = [(-t86 * t94 * t99 + t92 * t114) * t117, 0, 0, 0, 0, 0; (0.2e1 * (-t104 * t67 + t63 * t124) / t66 ^ 2 * (-t62 * t69 * t97 + t117 * t124) * t93 + ((0.2e1 * t102 * t63 * t69 - t123) * t62 + (-t63 * t123 + (-t67 + (-t92 * t131 - t116) * t84 * t123 - (t93 ^ 2 * t98 * t114 + (-t126 + (0.2e1 * t97 - t98) * t86) * t125) * t68 * t85) * t102) * qJD(1)) / t66) * t99, 0, 0, 0, 0, 0; 0.2e1 * (t78 * t127 + t77 * t80) * t130 + (-(-t82 * qJD(1) + t81 * qJD(2)) * t77 - 0.2e1 * t127 * t128 + (-t80 * t71 - (t113 * qJD(1) - t80 * qJD(2)) * t82 + t81 * t70) * t78) * t73, -0.2e1 * t130 + 0.2e1 * (-t73 * t129 + (t73 * t128 - t78 * t130) * t82) * t82, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:15
	% EndTime: 2019-10-10 12:42:16
	% DurationCPUTime: 1.04s
	% Computational Cost: add. (1479->91), mult. (4303->201), div. (668->14), fcn. (5516->11), ass. (0->91)
	t171 = sin(qJ(2));
	t172 = sin(qJ(1));
	t225 = cos(pkin(6));
	t195 = t172 * t225;
	t193 = t171 * t195;
	t174 = cos(qJ(2));
	t175 = cos(qJ(1));
	t209 = t175 * t174;
	t157 = -t193 + t209;
	t170 = sin(qJ(3));
	t173 = cos(qJ(3));
	t169 = sin(pkin(6));
	t213 = t169 * t172;
	t185 = -t157 * t170 + t173 * t213;
	t227 = t185 * qJD(3);
	t194 = t175 * t225;
	t192 = t174 * t194;
	t210 = t172 * t171;
	t153 = -t192 + t210;
	t212 = t169 * t174;
	t147 = atan2(-t153, -t212);
	t145 = sin(t147);
	t146 = cos(t147);
	t151 = t153 ^ 2;
	t165 = 0.1e1 / t169 ^ 2;
	t167 = 0.1e1 / t174 ^ 2;
	t150 = t151 * t165 * t167 + 0.1e1;
	t148 = 0.1e1 / t150;
	t164 = 0.1e1 / t169;
	t166 = 0.1e1 / t174;
	t199 = t153 * t164 * t166;
	t226 = (t146 * t199 - t145) * t148 + t145;
	t129 = -t145 * t153 - t146 * t212;
	t126 = 0.1e1 / t129;
	t144 = t157 * t173 + t170 * t213;
	t138 = 0.1e1 / t144;
	t127 = 0.1e1 / t129 ^ 2;
	t139 = 0.1e1 / t144 ^ 2;
	t182 = -t171 * t194 - t172 * t174;
	t183 = -t175 * t171 - t174 * t195;
	t135 = -t183 * qJD(1) - t182 * qJD(2);
	t207 = qJD(2) * t171;
	t196 = t167 * t207;
	t184 = t135 * t166 + t153 * t196;
	t215 = t148 * t164;
	t118 = t184 * t215;
	t188 = t145 * t212 - t146 * t153;
	t200 = t146 * t169 * t171;
	t114 = qJD(2) * t200 + t188 * t118 - t145 * t135;
	t224 = t114 * t126 * t127;
	t214 = t167 * t171;
	t187 = t153 * t214 - t166 * t182;
	t119 = t187 * t215;
	t115 = t188 * t119 + t145 * t182 + t200;
	t223 = t115 * t183;
	t134 = t182 * qJD(1) + t183 * qJD(2);
	t208 = qJD(1) * t169;
	t197 = t175 * t208;
	t124 = t144 * qJD(3) + t134 * t170 - t173 * t197;
	t137 = t185 ^ 2;
	t132 = t137 * t139 + 0.1e1;
	t218 = t139 * t185;
	t125 = t134 * t173 + t170 * t197 + t227;
	t220 = t125 * t138 * t139;
	t222 = (-t124 * t218 - t137 * t220) / t132 ^ 2;
	t168 = t166 * t167;
	t221 = (t135 * t153 * t167 + t151 * t168 * t207) * t165 / t150 ^ 2;
	t191 = qJD(2) * t225 + qJD(1);
	t206 = qJD(2) * t174;
	t133 = -qJD(1) * t192 - t175 * t206 + t191 * t210;
	t219 = t133 * t127;
	t217 = t145 * t183;
	t216 = t146 * t183;
	t211 = t169 * t175;
	t152 = t183 ^ 2;
	t122 = t152 * t127 + 0.1e1;
	t205 = 0.2e1 * (-t152 * t224 + t183 * t219) / t122 ^ 2;
	t204 = 0.2e1 * t224;
	t203 = 0.2e1 * t222;
	t202 = -0.2e1 * t221;
	t201 = t185 * t220;
	t198 = t172 * t208;
	t189 = t170 * t138 + t173 * t218;
	t186 = -t170 * t182 + t173 * t211;
	t142 = t170 * t211 + t173 * t182;
	t136 = -qJD(1) * t193 - t172 * t207 + t191 * t209;
	t130 = 0.1e1 / t132;
	t120 = 0.1e1 / t122;
	t117 = t226 * t183;
	t113 = (t187 * t202 + (t135 * t214 + t136 * t166 + (-t182 * t214 + (0.2e1 * t168 * t171 ^ 2 + t166) * t153) * qJD(2)) * t148) * t164;
	t1 = [(-t183 * t166 * t202 + (-t133 * t166 - t183 * t196) * t148) * t164, t113, 0, 0, 0, 0; t153 * t126 * t205 + (-t135 * t126 + (t114 * t153 + t117 * t133) * t127) * t120 - ((t117 * t204 - t226 * t219) * t120 + (t117 * t205 + ((t118 * t148 * t199 + t202) * t217 + (0.2e1 * t199 * t221 - t118 + (-t184 * t164 + t118) * t148) * t216) * t120) * t127) * t183, (-t126 * t157 - t127 * t223) * t205 + (-t204 * t223 + t134 * t126 + (-t157 * t114 + t115 * t133 + (t169 * t206 - t113 * t153 - t119 * t135 + (t119 * t212 + t182) * t118) * t216 + (t118 * t119 * t153 - t136 + (t113 * t174 + (-qJD(2) * t119 - t118) * t171) * t169) * t217) * t127) * t120, 0, 0, 0, 0; (t138 * t186 - t142 * t218) * t203 + ((t142 * qJD(3) - t136 * t170 + t173 * t198) * t138 - 0.2e1 * t142 * t201 + (t186 * t125 + (t186 * qJD(3) - t136 * t173 - t170 * t198) * t185 - t142 * t124) * t139) * t130, -t189 * t183 * t203 + (t189 * t133 - ((-qJD(3) * t138 + 0.2e1 * t201) * t173 + (t124 * t173 + (t125 + t227) * t170) * t139) * t183) * t130, -0.2e1 * t222 - 0.2e1 * (t124 * t139 * t130 - (-t130 * t220 - t139 * t222) * t185) * t185, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:15
	% EndTime: 2019-10-10 12:42:16
	% DurationCPUTime: 1.00s
	% Computational Cost: add. (2050->93), mult. (4665->202), div. (686->14), fcn. (5929->11), ass. (0->95)
	t202 = cos(qJ(2));
	t203 = cos(qJ(1));
	t253 = cos(pkin(6));
	t224 = t203 * t253;
	t222 = t202 * t224;
	t200 = sin(qJ(2));
	t201 = sin(qJ(1));
	t238 = t201 * t200;
	t179 = -t222 + t238;
	t199 = sin(pkin(6));
	t240 = t199 * t202;
	t173 = atan2(-t179, -t240);
	t171 = sin(t173);
	t172 = cos(t173);
	t177 = t179 ^ 2;
	t193 = 0.1e1 / t199 ^ 2;
	t196 = 0.1e1 / t202 ^ 2;
	t176 = t177 * t193 * t196 + 0.1e1;
	t174 = 0.1e1 / t176;
	t192 = 0.1e1 / t199;
	t195 = 0.1e1 / t202;
	t227 = t179 * t192 * t195;
	t254 = t174 * (t172 * t227 - t171) + t171;
	t155 = -t171 * t179 - t172 * t240;
	t152 = 0.1e1 / t155;
	t225 = t201 * t253;
	t223 = t200 * t225;
	t237 = t203 * t202;
	t183 = -t223 + t237;
	t198 = qJ(3) + qJ(4);
	t190 = sin(t198);
	t191 = cos(t198);
	t241 = t199 * t201;
	t170 = t183 * t191 + t190 * t241;
	t164 = 0.1e1 / t170;
	t153 = 0.1e1 / t155 ^ 2;
	t165 = 0.1e1 / t170 ^ 2;
	t210 = -t200 * t224 - t201 * t202;
	t211 = -t203 * t200 - t202 * t225;
	t161 = -qJD(1) * t211 - qJD(2) * t210;
	t235 = qJD(2) * t200;
	t226 = t196 * t235;
	t212 = t161 * t195 + t179 * t226;
	t243 = t174 * t192;
	t144 = t212 * t243;
	t216 = t171 * t240 - t172 * t179;
	t228 = t172 * t199 * t200;
	t140 = qJD(2) * t228 + t144 * t216 - t171 * t161;
	t252 = t140 * t152 * t153;
	t242 = t196 * t200;
	t215 = t179 * t242 - t195 * t210;
	t145 = t215 * t243;
	t141 = t145 * t216 + t171 * t210 + t228;
	t251 = t141 * t211;
	t194 = qJD(3) + qJD(4);
	t236 = qJD(1) * t199;
	t213 = -t183 * t194 + t203 * t236;
	t160 = qJD(1) * t210 + qJD(2) * t211;
	t221 = t194 * t241 + t160;
	t148 = t190 * t221 - t191 * t213;
	t169 = t183 * t190 - t191 * t241;
	t163 = t169 ^ 2;
	t158 = t163 * t165 + 0.1e1;
	t246 = t165 * t169;
	t149 = t190 * t213 + t191 * t221;
	t249 = t149 * t164 * t165;
	t250 = (t148 * t246 - t163 * t249) / t158 ^ 2;
	t197 = t195 * t196;
	t248 = (t161 * t179 * t196 + t177 * t197 * t235) * t193 / t176 ^ 2;
	t219 = qJD(2) * t253 + qJD(1);
	t234 = qJD(2) * t202;
	t159 = -qJD(1) * t222 - t203 * t234 + t219 * t238;
	t247 = t159 * t153;
	t245 = t171 * t211;
	t244 = t172 * t211;
	t239 = t199 * t203;
	t178 = t211 ^ 2;
	t150 = t178 * t153 + 0.1e1;
	t233 = 0.2e1 * (-t178 * t252 + t211 * t247) / t150 ^ 2;
	t232 = 0.2e1 * t252;
	t231 = 0.2e1 * t250;
	t230 = -0.2e1 * t248;
	t229 = t169 * t249;
	t162 = -qJD(1) * t223 - t201 * t235 + t219 * t237;
	t220 = t194 * t239 - t162;
	t217 = t190 * t164 - t191 * t246;
	t214 = t194 * t210 + t201 * t236;
	t168 = t190 * t239 + t191 * t210;
	t167 = t190 * t210 - t191 * t239;
	t156 = 0.1e1 / t158;
	t146 = 0.1e1 / t150;
	t143 = t254 * t211;
	t139 = (t215 * t230 + (t161 * t242 + t162 * t195 + (-t210 * t242 + (0.2e1 * t197 * t200 ^ 2 + t195) * t179) * qJD(2)) * t174) * t192;
	t137 = -0.2e1 * t250 + 0.2e1 * (t148 * t165 * t156 + (-t156 * t249 - t165 * t250) * t169) * t169;
	t1 = [(-t211 * t195 * t230 + (-t159 * t195 - t211 * t226) * t174) * t192, t139, 0, 0, 0, 0; t179 * t152 * t233 + (-t161 * t152 + (t140 * t179 + t143 * t159) * t153) * t146 - ((t143 * t232 - t254 * t247) * t146 + (t143 * t233 + ((t144 * t174 * t227 + t230) * t245 + (0.2e1 * t227 * t248 - t144 + (-t192 * t212 + t144) * t174) * t244) * t146) * t153) * t211, (-t152 * t183 - t153 * t251) * t233 + (-t232 * t251 + t160 * t152 + (-t183 * t140 + t141 * t159 + (t199 * t234 - t139 * t179 - t145 * t161 + (t145 * t240 + t210) * t144) * t244 + (t144 * t145 * t179 - t162 + (t139 * t202 + (-qJD(2) * t145 - t144) * t200) * t199) * t245) * t153) * t146, 0, 0, 0, 0; (-t164 * t167 + t168 * t246) * t231 + ((t190 * t220 + t191 * t214) * t164 + 0.2e1 * t168 * t229 + (-t167 * t149 - (-t190 * t214 + t191 * t220) * t169 - t168 * t148) * t165) * t156, -t217 * t211 * t231 + (t217 * t159 - ((-t164 * t194 - 0.2e1 * t229) * t191 + (t148 * t191 + (-t169 * t194 + t149) * t190) * t165) * t211) * t156, t137, t137, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:15
	% EndTime: 2019-10-10 12:42:16
	% DurationCPUTime: 0.95s
	% Computational Cost: add. (2443->93), mult. (4665->202), div. (686->14), fcn. (5929->11), ass. (0->95)
	t207 = cos(qJ(2));
	t208 = cos(qJ(1));
	t258 = cos(pkin(6));
	t229 = t208 * t258;
	t227 = t207 * t229;
	t205 = sin(qJ(2));
	t206 = sin(qJ(1));
	t243 = t206 * t205;
	t184 = -t227 + t243;
	t204 = sin(pkin(6));
	t245 = t204 * t207;
	t178 = atan2(-t184, -t245);
	t176 = sin(t178);
	t177 = cos(t178);
	t182 = t184 ^ 2;
	t199 = 0.1e1 / t204 ^ 2;
	t202 = 0.1e1 / t207 ^ 2;
	t181 = t182 * t199 * t202 + 0.1e1;
	t179 = 0.1e1 / t181;
	t198 = 0.1e1 / t204;
	t201 = 0.1e1 / t207;
	t232 = t184 * t198 * t201;
	t259 = t179 * (t177 * t232 - t176) + t176;
	t160 = -t176 * t184 - t177 * t245;
	t157 = 0.1e1 / t160;
	t230 = t206 * t258;
	t228 = t205 * t230;
	t242 = t208 * t207;
	t188 = -t228 + t242;
	t197 = qJ(3) + qJ(4) + pkin(12);
	t195 = sin(t197);
	t196 = cos(t197);
	t246 = t204 * t206;
	t175 = t188 * t196 + t195 * t246;
	t165 = 0.1e1 / t175;
	t158 = 0.1e1 / t160 ^ 2;
	t166 = 0.1e1 / t175 ^ 2;
	t215 = -t205 * t229 - t206 * t207;
	t216 = -t208 * t205 - t207 * t230;
	t170 = -qJD(1) * t216 - qJD(2) * t215;
	t240 = qJD(2) * t205;
	t231 = t202 * t240;
	t217 = t170 * t201 + t184 * t231;
	t248 = t179 * t198;
	t149 = t217 * t248;
	t221 = t176 * t245 - t177 * t184;
	t233 = t177 * t204 * t205;
	t145 = qJD(2) * t233 + t149 * t221 - t176 * t170;
	t257 = t145 * t157 * t158;
	t200 = qJD(3) + qJD(4);
	t241 = qJD(1) * t204;
	t218 = -t188 * t200 + t208 * t241;
	t169 = qJD(1) * t215 + qJD(2) * t216;
	t226 = t200 * t246 + t169;
	t150 = t195 * t226 - t196 * t218;
	t174 = t188 * t195 - t196 * t246;
	t164 = t174 ^ 2;
	t163 = t164 * t166 + 0.1e1;
	t252 = t166 * t174;
	t151 = t195 * t218 + t196 * t226;
	t254 = t151 * t165 * t166;
	t256 = (t150 * t252 - t164 * t254) / t163 ^ 2;
	t247 = t202 * t205;
	t220 = t184 * t247 - t201 * t215;
	t152 = t220 * t248;
	t147 = t152 * t221 + t176 * t215 + t233;
	t255 = t147 * t216;
	t203 = t201 * t202;
	t253 = (t170 * t184 * t202 + t182 * t203 * t240) * t199 / t181 ^ 2;
	t224 = qJD(2) * t258 + qJD(1);
	t239 = qJD(2) * t207;
	t168 = -qJD(1) * t227 - t208 * t239 + t224 * t243;
	t251 = t168 * t158;
	t250 = t176 * t216;
	t249 = t177 * t216;
	t244 = t204 * t208;
	t183 = t216 ^ 2;
	t155 = t183 * t158 + 0.1e1;
	t238 = 0.2e1 * (-t183 * t257 + t216 * t251) / t155 ^ 2;
	t237 = 0.2e1 * t257;
	t236 = 0.2e1 * t256;
	t235 = -0.2e1 * t253;
	t234 = t174 * t254;
	t171 = -qJD(1) * t228 - t206 * t240 + t224 * t242;
	t225 = t200 * t244 - t171;
	t222 = t195 * t165 - t196 * t252;
	t219 = t200 * t215 + t206 * t241;
	t173 = t195 * t244 + t196 * t215;
	t172 = t195 * t215 - t196 * t244;
	t161 = 0.1e1 / t163;
	t153 = 0.1e1 / t155;
	t148 = t259 * t216;
	t144 = (t220 * t235 + (t170 * t247 + t171 * t201 + (-t215 * t247 + (0.2e1 * t203 * t205 ^ 2 + t201) * t184) * qJD(2)) * t179) * t198;
	t142 = -0.2e1 * t256 + 0.2e1 * (t150 * t166 * t161 + (-t161 * t254 - t166 * t256) * t174) * t174;
	t1 = [(-t216 * t201 * t235 + (-t168 * t201 - t216 * t231) * t179) * t198, t144, 0, 0, 0, 0; t184 * t157 * t238 + (-t170 * t157 + (t145 * t184 + t148 * t168) * t158) * t153 - ((t148 * t237 - t259 * t251) * t153 + (t148 * t238 + ((t149 * t179 * t232 + t235) * t250 + (0.2e1 * t232 * t253 - t149 + (-t198 * t217 + t149) * t179) * t249) * t153) * t158) * t216, (-t157 * t188 - t158 * t255) * t238 + (-t237 * t255 + t169 * t157 + (-t188 * t145 + t147 * t168 + (t204 * t239 - t144 * t184 - t152 * t170 + (t152 * t245 + t215) * t149) * t249 + (t149 * t152 * t184 - t171 + (t144 * t207 + (-qJD(2) * t152 - t149) * t205) * t204) * t250) * t158) * t153, 0, 0, 0, 0; (-t165 * t172 + t173 * t252) * t236 + ((t195 * t225 + t196 * t219) * t165 + 0.2e1 * t173 * t234 + (-t172 * t151 - (-t195 * t219 + t196 * t225) * t174 - t173 * t150) * t166) * t161, -t222 * t216 * t236 + (t222 * t168 - ((-t165 * t200 - 0.2e1 * t234) * t196 + (t150 * t196 + (-t174 * t200 + t151) * t195) * t166) * t216) * t161, t142, t142, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:15
	% EndTime: 2019-10-10 12:42:18
	% DurationCPUTime: 2.56s
	% Computational Cost: add. (17252->153), mult. (18084->300), div. (989->12), fcn. (22873->13), ass. (0->136)
	t305 = qJ(3) + qJ(4) + pkin(12);
	t303 = sin(t305);
	t308 = cos(pkin(6));
	t312 = cos(qJ(2));
	t384 = sin(qJ(1));
	t341 = t384 * t312;
	t310 = sin(qJ(2));
	t313 = cos(qJ(1));
	t361 = t313 * t310;
	t324 = -t308 * t361 - t341;
	t304 = cos(t305);
	t307 = sin(pkin(6));
	t364 = t307 * t313;
	t345 = t304 * t364;
	t276 = -t303 * t324 + t345;
	t366 = t307 * t310;
	t347 = t303 * t366;
	t286 = -t308 * t304 + t347;
	t257 = atan2(-t276, t286);
	t252 = sin(t257);
	t253 = cos(t257);
	t243 = -t252 * t276 + t253 * t286;
	t241 = 0.1e1 / t243 ^ 2;
	t342 = t384 * t310;
	t334 = t308 * t342;
	t360 = t313 * t312;
	t294 = -t334 + t360;
	t343 = t307 * t384;
	t281 = t294 * t303 - t304 * t343;
	t271 = t281 ^ 2;
	t237 = t271 * t241 + 0.1e1;
	t323 = -t308 * t341 - t361;
	t273 = qJD(1) * t324 + qJD(2) * t323;
	t306 = qJD(3) + qJD(4);
	t330 = t306 * t343 + t273;
	t340 = qJD(1) * t364;
	t367 = t304 * t306;
	t247 = t294 * t367 + t303 * t330 - t304 * t340;
	t377 = t247 * t241;
	t270 = t276 ^ 2;
	t284 = 0.1e1 / t286 ^ 2;
	t256 = t270 * t284 + 0.1e1;
	t254 = 0.1e1 / t256;
	t338 = t384 * qJD(2);
	t275 = -qJD(1) * t334 - t310 * t338 + (qJD(2) * t308 + qJD(1)) * t360;
	t298 = t303 * t364;
	t339 = t384 * qJD(1);
	t333 = t307 * t339;
	t249 = t275 * t303 - t306 * t298 - t304 * t333 - t324 * t367;
	t358 = qJD(2) * t312;
	t326 = t306 * t308 + t307 * t358;
	t346 = t304 * t366;
	t268 = t303 * t326 + t306 * t346;
	t283 = 0.1e1 / t286;
	t371 = t276 * t284;
	t329 = -t249 * t283 + t268 * t371;
	t231 = t329 * t254;
	t331 = -t252 * t286 - t253 * t276;
	t226 = t231 * t331 - t252 * t249 + t253 * t268;
	t240 = 0.1e1 / t243;
	t242 = t240 * t241;
	t382 = t226 * t242;
	t356 = 0.2e1 * (-t271 * t382 + t281 * t377) / t237 ^ 2;
	t388 = t268 * t284;
	t344 = t308 * t360;
	t291 = -t342 + t344;
	t365 = t307 * t312;
	t325 = -t283 * t291 + t365 * t371;
	t387 = t303 * t325;
	t250 = t303 * (t306 * t324 + t333) + t275 * t304 - t306 * t345;
	t282 = t294 * t304 + t303 * t343;
	t311 = cos(qJ(6));
	t309 = sin(qJ(6));
	t369 = t323 * t309;
	t265 = t282 * t311 - t369;
	t259 = 0.1e1 / t265;
	t260 = 0.1e1 / t265 ^ 2;
	t386 = -0.2e1 * t276;
	t385 = 0.2e1 * t281;
	t248 = t330 * t304 + (-t294 * t306 + t340) * t303;
	t272 = -qJD(1) * t344 - t313 * t358 + (t308 * t338 + t339) * t310;
	t238 = qJD(6) * t265 + t248 * t309 + t272 * t311;
	t368 = t323 * t311;
	t264 = t282 * t309 + t368;
	t258 = t264 ^ 2;
	t246 = t258 * t260 + 0.1e1;
	t374 = t260 * t264;
	t357 = qJD(6) * t264;
	t239 = t248 * t311 - t272 * t309 - t357;
	t379 = t239 * t259 * t260;
	t381 = (t238 * t374 - t258 * t379) / t246 ^ 2;
	t373 = t283 * t388;
	t380 = (t249 * t371 - t270 * t373) / t256 ^ 2;
	t378 = t241 * t281;
	t376 = t252 * t281;
	t375 = t253 * t281;
	t372 = t276 * t283;
	t370 = t323 * t303;
	t363 = t309 * t259;
	t362 = t311 * t264;
	t359 = qJD(2) * t310;
	t355 = -0.2e1 * t381;
	t354 = 0.2e1 * t381;
	t353 = -0.2e1 * t380;
	t352 = t242 * t385;
	t351 = t283 * t380;
	t350 = t241 * t376;
	t349 = t241 * t375;
	t348 = t264 * t379;
	t337 = 0.2e1 * t348;
	t336 = t373 * t386;
	t278 = -t304 * t324 - t298;
	t332 = -qJD(6) * t304 * t323 + t273;
	t263 = -t278 * t311 + t291 * t309;
	t262 = -t278 * t309 - t291 * t311;
	t328 = t260 * t362 - t363;
	t287 = t308 * t303 + t346;
	t327 = -t278 * t283 + t287 * t371;
	t321 = -t252 + (t253 * t372 + t252) * t254;
	t320 = qJD(6) * t294 + t272 * t304 - t306 * t370;
	t274 = qJD(1) * t323 + qJD(2) * t324;
	t269 = t304 * t326 - t306 * t347;
	t267 = t294 * t309 + t304 * t368;
	t266 = -t294 * t311 + t304 * t369;
	t244 = 0.1e1 / t246;
	t235 = 0.1e1 / t237;
	t234 = t254 * t387;
	t232 = t327 * t254;
	t230 = t321 * t281;
	t228 = (-t252 * t291 + t253 * t365) * t303 + t331 * t234;
	t227 = t232 * t331 - t252 * t278 + t253 * t287;
	t224 = t327 * t353 + (t287 * t336 - t250 * t283 + (t249 * t287 + t268 * t278 + t269 * t276) * t284) * t254;
	t223 = t353 * t387 + (t325 * t367 + (t336 * t365 - t274 * t283 + (t268 * t291 + (t249 * t312 - t276 * t359) * t307) * t284) * t303) * t254;
	t222 = t328 * t281 * t355 + (t328 * t247 + ((-qJD(6) * t259 - 0.2e1 * t348) * t311 + (t238 * t311 + (t239 - t357) * t309) * t260) * t281) * t244;
	t221 = (t227 * t378 - t240 * t282) * t356 + (t227 * t226 * t352 + t248 * t240 + (-t282 * t226 - t227 * t247 - (-t224 * t276 - t232 * t249 + t269 + (-t232 * t286 - t278) * t231) * t375 - (-t224 * t286 - t232 * t268 - t250 + (t232 * t276 - t287) * t231) * t376) * t241) * t235;
	t1 = [t351 * t385 + (-t247 * t283 + t281 * t388) * t254, t223, t224, t224, 0, 0; t276 * t240 * t356 + (-t249 * t240 + (t226 * t276 - t230 * t247) * t241) * t235 + (t230 * t241 * t356 + (0.2e1 * t230 * t382 - (-t231 * t254 * t372 + t353) * t350 - (t351 * t386 - t231 + (t231 - t329) * t254) * t349 - t321 * t377) * t235) * t281, (t228 * t378 - t240 * t370) * t356 + (-t228 * t377 + (t272 * t303 + t323 * t367) * t240 + (t228 * t352 - t241 * t370) * t226 - (-t223 * t276 - t234 * t249 + (-t303 * t359 + t312 * t367) * t307 + (-t234 * t286 - t291 * t303) * t231) * t349 - (-t291 * t367 - t223 * t286 - t234 * t268 - t274 * t303 + (t234 * t276 - t303 * t365) * t231) * t350) * t235, t221, t221, 0, 0; (-t259 * t262 + t263 * t374) * t354 + ((qJD(6) * t263 - t250 * t309 - t274 * t311) * t259 + t263 * t337 + (-t262 * t239 - (-qJD(6) * t262 - t250 * t311 + t274 * t309) * t264 - t263 * t238) * t260) * t244, (-t259 * t266 + t267 * t374) * t354 + (t267 * t337 - t332 * t259 * t311 + t320 * t363 + (-t264 * t309 * t332 - t267 * t238 - t266 * t239 - t320 * t362) * t260) * t244, t222, t222, 0, t355 + 0.2e1 * (t238 * t260 * t244 + (-t244 * t379 - t260 * t381) * t264) * t264;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end