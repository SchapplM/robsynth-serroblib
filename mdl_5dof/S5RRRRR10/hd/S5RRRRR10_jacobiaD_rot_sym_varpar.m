% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRR10
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
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRRRR10_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 21:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRRR10_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR10_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:01:26
	% EndTime: 2019-12-29 21:01:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:01:26
	% EndTime: 2019-12-29 21:01:26
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:01:26
	% EndTime: 2019-12-29 21:01:27
	% DurationCPUTime: 0.62s
	% Computational Cost: add. (215->39), mult. (853->106), div. (126->12), fcn. (1047->9), ass. (0->54)
	t99 = sin(pkin(5));
	t93 = t99 ^ 2;
	t100 = cos(pkin(5));
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
	t1 = [(-t86 * t94 * t99 + t92 * t114) * t117, 0, 0, 0, 0; (0.2e1 * (-t104 * t67 + t63 * t124) / t66 ^ 2 * (-t62 * t69 * t97 + t117 * t124) * t93 + ((0.2e1 * t102 * t63 * t69 - t123) * t62 + (-t63 * t123 + (-t67 + (-t92 * t131 - t116) * t84 * t123 - (t93 ^ 2 * t98 * t114 + (-t126 + (0.2e1 * t97 - t98) * t86) * t125) * t68 * t85) * t102) * qJD(1)) / t66) * t99, 0, 0, 0, 0; 0.2e1 * (t78 * t127 + t77 * t80) * t130 + (-(-t82 * qJD(1) + t81 * qJD(2)) * t77 - 0.2e1 * t127 * t128 + (-t80 * t71 - (t113 * qJD(1) - t80 * qJD(2)) * t82 + t81 * t70) * t78) * t73, -0.2e1 * t130 + 0.2e1 * (-t73 * t129 + (t73 * t128 - t78 * t130) * t82) * t82, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:01:26
	% EndTime: 2019-12-29 21:01:28
	% DurationCPUTime: 1.58s
	% Computational Cost: add. (1479->91), mult. (4303->201), div. (668->14), fcn. (5516->11), ass. (0->91)
	t171 = sin(qJ(2));
	t172 = sin(qJ(1));
	t225 = cos(pkin(5));
	t195 = t172 * t225;
	t193 = t171 * t195;
	t174 = cos(qJ(2));
	t175 = cos(qJ(1));
	t209 = t175 * t174;
	t157 = -t193 + t209;
	t170 = sin(qJ(3));
	t173 = cos(qJ(3));
	t169 = sin(pkin(5));
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
	t122 = t127 * t152 + 0.1e1;
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
	t1 = [(-t183 * t166 * t202 + (-t133 * t166 - t183 * t196) * t148) * t164, t113, 0, 0, 0; t153 * t126 * t205 + (-t135 * t126 + (t114 * t153 + t117 * t133) * t127) * t120 - ((t117 * t204 - t226 * t219) * t120 + (t117 * t205 + ((t118 * t148 * t199 + t202) * t217 + (0.2e1 * t199 * t221 - t118 + (-t184 * t164 + t118) * t148) * t216) * t120) * t127) * t183, (-t126 * t157 - t127 * t223) * t205 + (-t204 * t223 + t134 * t126 + (-t157 * t114 + t115 * t133 + (t169 * t206 - t113 * t153 - t119 * t135 + (t119 * t212 + t182) * t118) * t216 + (t118 * t119 * t153 - t136 + (t113 * t174 + (-qJD(2) * t119 - t118) * t171) * t169) * t217) * t127) * t120, 0, 0, 0; (t138 * t186 - t142 * t218) * t203 + ((t142 * qJD(3) - t136 * t170 + t173 * t198) * t138 - 0.2e1 * t142 * t201 + (t186 * t125 + (t186 * qJD(3) - t136 * t173 - t170 * t198) * t185 - t142 * t124) * t139) * t130, -t189 * t183 * t203 + (t189 * t133 - ((-qJD(3) * t138 + 0.2e1 * t201) * t173 + (t124 * t173 + (t125 + t227) * t170) * t139) * t183) * t130, -0.2e1 * t222 - 0.2e1 * (t124 * t139 * t130 - (-t130 * t220 - t139 * t222) * t185) * t185, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:01:27
	% EndTime: 2019-12-29 21:01:28
	% DurationCPUTime: 1.45s
	% Computational Cost: add. (2050->93), mult. (4665->202), div. (686->14), fcn. (5929->11), ass. (0->95)
	t202 = cos(qJ(2));
	t203 = cos(qJ(1));
	t253 = cos(pkin(5));
	t224 = t203 * t253;
	t222 = t202 * t224;
	t200 = sin(qJ(2));
	t201 = sin(qJ(1));
	t238 = t201 * t200;
	t179 = -t222 + t238;
	t199 = sin(pkin(5));
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
	t1 = [(-t211 * t195 * t230 + (-t159 * t195 - t211 * t226) * t174) * t192, t139, 0, 0, 0; t179 * t152 * t233 + (-t161 * t152 + (t140 * t179 + t143 * t159) * t153) * t146 - ((t143 * t232 - t254 * t247) * t146 + (t143 * t233 + ((t144 * t174 * t227 + t230) * t245 + (0.2e1 * t227 * t248 - t144 + (-t192 * t212 + t144) * t174) * t244) * t146) * t153) * t211, (-t152 * t183 - t153 * t251) * t233 + (-t232 * t251 + t160 * t152 + (-t183 * t140 + t141 * t159 + (t199 * t234 - t139 * t179 - t145 * t161 + (t145 * t240 + t210) * t144) * t244 + (t144 * t145 * t179 - t162 + (t139 * t202 + (-qJD(2) * t145 - t144) * t200) * t199) * t245) * t153) * t146, 0, 0, 0; (-t164 * t167 + t168 * t246) * t231 + ((t190 * t220 + t191 * t214) * t164 + 0.2e1 * t168 * t229 + (-t167 * t149 - (-t190 * t214 + t191 * t220) * t169 - t168 * t148) * t165) * t156, -t217 * t211 * t231 + (t217 * t159 - ((-t164 * t194 - 0.2e1 * t229) * t191 + (t148 * t191 + (-t169 * t194 + t149) * t190) * t165) * t211) * t156, t137, t137, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:01:27
	% EndTime: 2019-12-29 21:01:31
	% DurationCPUTime: 3.84s
	% Computational Cost: add. (11947->153), mult. (18084->300), div. (989->12), fcn. (22873->13), ass. (0->136)
	t301 = qJ(3) + qJ(4);
	t298 = sin(t301);
	t303 = cos(pkin(5));
	t307 = cos(qJ(2));
	t379 = sin(qJ(1));
	t336 = t379 * t307;
	t305 = sin(qJ(2));
	t308 = cos(qJ(1));
	t356 = t308 * t305;
	t319 = -t303 * t356 - t336;
	t299 = cos(t301);
	t302 = sin(pkin(5));
	t358 = t302 * t308;
	t340 = t299 * t358;
	t271 = -t298 * t319 + t340;
	t360 = t302 * t305;
	t342 = t298 * t360;
	t281 = -t303 * t299 + t342;
	t260 = atan2(-t271, t281);
	t251 = sin(t260);
	t252 = cos(t260);
	t238 = -t251 * t271 + t252 * t281;
	t236 = 0.1e1 / t238 ^ 2;
	t337 = t379 * t305;
	t329 = t303 * t337;
	t355 = t308 * t307;
	t287 = -t329 + t355;
	t338 = t302 * t379;
	t276 = t287 * t298 - t299 * t338;
	t270 = t276 ^ 2;
	t232 = t270 * t236 + 0.1e1;
	t318 = -t303 * t336 - t356;
	t266 = t319 * qJD(1) + t318 * qJD(2);
	t300 = qJD(3) + qJD(4);
	t325 = t300 * t338 + t266;
	t335 = qJD(1) * t358;
	t361 = t299 * t300;
	t242 = t287 * t361 + t325 * t298 - t299 * t335;
	t372 = t236 * t276;
	t269 = t271 ^ 2;
	t279 = 0.1e1 / t281 ^ 2;
	t259 = t269 * t279 + 0.1e1;
	t257 = 0.1e1 / t259;
	t334 = qJD(2) * t379;
	t268 = -qJD(1) * t329 - t305 * t334 + (qJD(2) * t303 + qJD(1)) * t355;
	t293 = t298 * t358;
	t333 = t379 * qJD(1);
	t328 = t302 * t333;
	t244 = t268 * t298 - t300 * t293 - t299 * t328 - t319 * t361;
	t353 = qJD(2) * t307;
	t321 = t300 * t303 + t302 * t353;
	t341 = t299 * t360;
	t263 = t321 * t298 + t300 * t341;
	t278 = 0.1e1 / t281;
	t365 = t271 * t279;
	t324 = -t244 * t278 + t263 * t365;
	t226 = t324 * t257;
	t326 = -t251 * t281 - t252 * t271;
	t221 = t326 * t226 - t251 * t244 + t252 * t263;
	t235 = 0.1e1 / t238;
	t237 = t235 * t236;
	t377 = t221 * t237;
	t351 = 0.2e1 * (t242 * t372 - t270 * t377) / t232 ^ 2;
	t383 = t263 * t279;
	t339 = t303 * t355;
	t284 = -t337 + t339;
	t359 = t302 * t307;
	t320 = -t278 * t284 + t359 * t365;
	t382 = t298 * t320;
	t245 = (t300 * t319 + t328) * t298 + t268 * t299 - t300 * t340;
	t277 = t287 * t299 + t298 * t338;
	t306 = cos(qJ(5));
	t304 = sin(qJ(5));
	t363 = t318 * t304;
	t256 = t277 * t306 - t363;
	t248 = 0.1e1 / t256;
	t249 = 0.1e1 / t256 ^ 2;
	t381 = -0.2e1 * t271;
	t380 = 0.2e1 * t276;
	t243 = t325 * t299 + (-t287 * t300 + t335) * t298;
	t265 = -qJD(1) * t339 - t308 * t353 + (t303 * t334 + t333) * t305;
	t233 = t256 * qJD(5) + t243 * t304 + t265 * t306;
	t362 = t318 * t306;
	t255 = t277 * t304 + t362;
	t247 = t255 ^ 2;
	t241 = t247 * t249 + 0.1e1;
	t371 = t249 * t255;
	t352 = qJD(5) * t255;
	t234 = t243 * t306 - t265 * t304 - t352;
	t374 = t234 * t248 * t249;
	t376 = (t233 * t371 - t247 * t374) / t241 ^ 2;
	t367 = t278 * t383;
	t375 = (t244 * t365 - t269 * t367) / t259 ^ 2;
	t373 = t236 * t242;
	t370 = t251 * t276;
	t369 = t252 * t276;
	t368 = t255 * t306;
	t366 = t271 * t278;
	t364 = t318 * t298;
	t357 = t304 * t248;
	t354 = qJD(2) * t305;
	t350 = -0.2e1 * t376;
	t349 = 0.2e1 * t376;
	t348 = -0.2e1 * t375;
	t347 = t237 * t380;
	t346 = t278 * t375;
	t345 = t236 * t370;
	t344 = t236 * t369;
	t343 = t255 * t374;
	t332 = 0.2e1 * t343;
	t331 = t367 * t381;
	t273 = -t299 * t319 - t293;
	t327 = -qJD(5) * t299 * t318 + t266;
	t254 = -t273 * t306 + t284 * t304;
	t253 = -t273 * t304 - t284 * t306;
	t323 = t249 * t368 - t357;
	t282 = t303 * t298 + t341;
	t322 = -t273 * t278 + t282 * t365;
	t316 = -t251 + (t252 * t366 + t251) * t257;
	t315 = qJD(5) * t287 + t265 * t299 - t300 * t364;
	t267 = t318 * qJD(1) + t319 * qJD(2);
	t264 = t321 * t299 - t300 * t342;
	t262 = t287 * t304 + t299 * t362;
	t261 = -t287 * t306 + t299 * t363;
	t239 = 0.1e1 / t241;
	t230 = 0.1e1 / t232;
	t229 = t257 * t382;
	t228 = t322 * t257;
	t225 = t316 * t276;
	t223 = (-t251 * t284 + t252 * t359) * t298 + t326 * t229;
	t222 = t326 * t228 - t251 * t273 + t252 * t282;
	t219 = t322 * t348 + (t282 * t331 - t245 * t278 + (t244 * t282 + t263 * t273 + t264 * t271) * t279) * t257;
	t218 = t348 * t382 + (t320 * t361 + (t331 * t359 - t267 * t278 + (t263 * t284 + (t244 * t307 - t271 * t354) * t302) * t279) * t298) * t257;
	t217 = t323 * t276 * t350 + (t323 * t242 + ((-qJD(5) * t248 - 0.2e1 * t343) * t306 + (t233 * t306 + (t234 - t352) * t304) * t249) * t276) * t239;
	t216 = (t222 * t372 - t235 * t277) * t351 + (t222 * t221 * t347 + t243 * t235 + (-t277 * t221 - t222 * t242 - (-t219 * t271 - t228 * t244 + t264 + (-t228 * t281 - t273) * t226) * t369 - (-t219 * t281 - t228 * t263 - t245 + (t228 * t271 - t282) * t226) * t370) * t236) * t230;
	t1 = [t346 * t380 + (-t242 * t278 + t276 * t383) * t257, t218, t219, t219, 0; t271 * t235 * t351 + (-t244 * t235 + (t221 * t271 - t225 * t242) * t236) * t230 + (t225 * t236 * t351 + (0.2e1 * t225 * t377 - (-t226 * t257 * t366 + t348) * t345 - (t346 * t381 - t226 + (t226 - t324) * t257) * t344 - t316 * t373) * t230) * t276, (t223 * t372 - t235 * t364) * t351 + (-t223 * t373 + (t265 * t298 + t318 * t361) * t235 + (t223 * t347 - t236 * t364) * t221 - (-t218 * t271 - t229 * t244 + (-t298 * t354 + t307 * t361) * t302 + (-t229 * t281 - t284 * t298) * t226) * t344 - (-t284 * t361 - t218 * t281 - t229 * t263 - t267 * t298 + (t229 * t271 - t298 * t359) * t226) * t345) * t230, t216, t216, 0; (-t248 * t253 + t254 * t371) * t349 + ((qJD(5) * t254 - t245 * t304 - t267 * t306) * t248 + t254 * t332 + (-t253 * t234 - (-qJD(5) * t253 - t245 * t306 + t267 * t304) * t255 - t254 * t233) * t249) * t239, (-t248 * t261 + t262 * t371) * t349 + (t262 * t332 - t327 * t248 * t306 + t315 * t357 + (-t255 * t304 * t327 - t262 * t233 - t261 * t234 - t315 * t368) * t249) * t239, t217, t217, t350 + 0.2e1 * (t233 * t249 * t239 + (-t239 * t374 - t249 * t376) * t255) * t255;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end