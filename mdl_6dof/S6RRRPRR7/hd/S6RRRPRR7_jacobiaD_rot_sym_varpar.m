% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR7
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
%   Wie in S6RRRPRR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:04
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
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
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:13
	% DurationCPUTime: 1.06s
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
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:13
	% DurationCPUTime: 1.04s
	% Computational Cost: add. (1788->92), mult. (4303->201), div. (668->14), fcn. (5516->11), ass. (0->92)
	t179 = sin(qJ(2));
	t180 = sin(qJ(1));
	t232 = cos(pkin(6));
	t202 = t180 * t232;
	t200 = t179 * t202;
	t181 = cos(qJ(2));
	t182 = cos(qJ(1));
	t216 = t182 * t181;
	t163 = -t200 + t216;
	t174 = qJ(3) + pkin(12);
	t170 = sin(t174);
	t171 = cos(t174);
	t178 = sin(pkin(6));
	t220 = t178 * t180;
	t192 = -t163 * t170 + t171 * t220;
	t234 = qJD(3) * t192;
	t201 = t182 * t232;
	t199 = t181 * t201;
	t217 = t180 * t179;
	t159 = -t199 + t217;
	t219 = t178 * t181;
	t153 = atan2(-t159, -t219);
	t151 = sin(t153);
	t152 = cos(t153);
	t157 = t159 ^ 2;
	t173 = 0.1e1 / t178 ^ 2;
	t176 = 0.1e1 / t181 ^ 2;
	t156 = t157 * t173 * t176 + 0.1e1;
	t154 = 0.1e1 / t156;
	t172 = 0.1e1 / t178;
	t175 = 0.1e1 / t181;
	t206 = t159 * t172 * t175;
	t233 = (t152 * t206 - t151) * t154 + t151;
	t135 = -t151 * t159 - t152 * t219;
	t132 = 0.1e1 / t135;
	t150 = t163 * t171 + t170 * t220;
	t144 = 0.1e1 / t150;
	t133 = 0.1e1 / t135 ^ 2;
	t145 = 0.1e1 / t150 ^ 2;
	t189 = -t179 * t201 - t180 * t181;
	t190 = -t182 * t179 - t181 * t202;
	t141 = -qJD(1) * t190 - qJD(2) * t189;
	t214 = qJD(2) * t179;
	t203 = t176 * t214;
	t191 = t141 * t175 + t159 * t203;
	t222 = t154 * t172;
	t124 = t191 * t222;
	t195 = t151 * t219 - t152 * t159;
	t207 = t152 * t178 * t179;
	t120 = qJD(2) * t207 + t124 * t195 - t151 * t141;
	t231 = t120 * t132 * t133;
	t221 = t176 * t179;
	t194 = t159 * t221 - t175 * t189;
	t125 = t194 * t222;
	t121 = t125 * t195 + t151 * t189 + t207;
	t230 = t121 * t190;
	t140 = qJD(1) * t189 + qJD(2) * t190;
	t215 = qJD(1) * t178;
	t204 = t182 * t215;
	t129 = qJD(3) * t150 + t140 * t170 - t171 * t204;
	t143 = t192 ^ 2;
	t138 = t143 * t145 + 0.1e1;
	t225 = t145 * t192;
	t130 = t140 * t171 + t170 * t204 + t234;
	t228 = t130 * t144 * t145;
	t229 = (-t129 * t225 - t143 * t228) / t138 ^ 2;
	t177 = t175 * t176;
	t227 = (t141 * t159 * t176 + t157 * t177 * t214) * t173 / t156 ^ 2;
	t198 = qJD(2) * t232 + qJD(1);
	t213 = qJD(2) * t181;
	t139 = -qJD(1) * t199 - t182 * t213 + t198 * t217;
	t226 = t139 * t133;
	t224 = t151 * t190;
	t223 = t152 * t190;
	t218 = t178 * t182;
	t158 = t190 ^ 2;
	t128 = t158 * t133 + 0.1e1;
	t212 = 0.2e1 * (-t158 * t231 + t190 * t226) / t128 ^ 2;
	t211 = 0.2e1 * t231;
	t210 = 0.2e1 * t229;
	t209 = -0.2e1 * t227;
	t208 = t192 * t228;
	t205 = t180 * t215;
	t196 = t170 * t144 + t171 * t225;
	t193 = -t170 * t189 + t171 * t218;
	t148 = t170 * t218 + t171 * t189;
	t142 = -qJD(1) * t200 - t180 * t214 + t198 * t216;
	t136 = 0.1e1 / t138;
	t126 = 0.1e1 / t128;
	t123 = t233 * t190;
	t119 = (t194 * t209 + (t141 * t221 + t142 * t175 + (-t189 * t221 + (0.2e1 * t177 * t179 ^ 2 + t175) * t159) * qJD(2)) * t154) * t172;
	t1 = [(-t190 * t175 * t209 + (-t139 * t175 - t190 * t203) * t154) * t172, t119, 0, 0, 0, 0; t159 * t132 * t212 + (-t141 * t132 + (t120 * t159 + t123 * t139) * t133) * t126 - ((t123 * t211 - t233 * t226) * t126 + (t123 * t212 + ((t124 * t154 * t206 + t209) * t224 + (0.2e1 * t206 * t227 - t124 + (-t172 * t191 + t124) * t154) * t223) * t126) * t133) * t190, (-t132 * t163 - t133 * t230) * t212 + (-t211 * t230 + t140 * t132 + (-t163 * t120 + t121 * t139 + (t178 * t213 - t119 * t159 - t125 * t141 + (t125 * t219 + t189) * t124) * t223 + (t124 * t125 * t159 - t142 + (t119 * t181 + (-qJD(2) * t125 - t124) * t179) * t178) * t224) * t133) * t126, 0, 0, 0, 0; (t144 * t193 - t148 * t225) * t210 + ((qJD(3) * t148 - t142 * t170 + t171 * t205) * t144 - 0.2e1 * t148 * t208 + (t193 * t130 + (qJD(3) * t193 - t142 * t171 - t170 * t205) * t192 - t148 * t129) * t145) * t136, -t196 * t190 * t210 + (t196 * t139 - ((-qJD(3) * t144 + 0.2e1 * t208) * t171 + (t129 * t171 + (t130 + t234) * t170) * t145) * t190) * t136, -0.2e1 * t229 - 0.2e1 * (t129 * t145 * t136 - (-t136 * t228 - t145 * t229) * t192) * t192, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:13
	% DurationCPUTime: 1.00s
	% Computational Cost: add. (2443->93), mult. (4665->202), div. (686->14), fcn. (5929->11), ass. (0->95)
	t205 = cos(qJ(2));
	t206 = cos(qJ(1));
	t256 = cos(pkin(6));
	t227 = t206 * t256;
	t225 = t205 * t227;
	t203 = sin(qJ(2));
	t204 = sin(qJ(1));
	t241 = t204 * t203;
	t182 = -t225 + t241;
	t202 = sin(pkin(6));
	t243 = t202 * t205;
	t176 = atan2(-t182, -t243);
	t174 = sin(t176);
	t175 = cos(t176);
	t180 = t182 ^ 2;
	t197 = 0.1e1 / t202 ^ 2;
	t200 = 0.1e1 / t205 ^ 2;
	t179 = t180 * t197 * t200 + 0.1e1;
	t177 = 0.1e1 / t179;
	t196 = 0.1e1 / t202;
	t199 = 0.1e1 / t205;
	t230 = t182 * t196 * t199;
	t257 = (t175 * t230 - t174) * t177 + t174;
	t158 = -t174 * t182 - t175 * t243;
	t155 = 0.1e1 / t158;
	t228 = t204 * t256;
	t226 = t203 * t228;
	t240 = t206 * t205;
	t186 = -t226 + t240;
	t195 = qJ(3) + pkin(12) + qJ(5);
	t193 = sin(t195);
	t194 = cos(t195);
	t244 = t202 * t204;
	t173 = t186 * t194 + t193 * t244;
	t163 = 0.1e1 / t173;
	t156 = 0.1e1 / t158 ^ 2;
	t164 = 0.1e1 / t173 ^ 2;
	t213 = -t203 * t227 - t204 * t205;
	t214 = -t206 * t203 - t205 * t228;
	t168 = -t214 * qJD(1) - t213 * qJD(2);
	t238 = qJD(2) * t203;
	t229 = t200 * t238;
	t215 = t168 * t199 + t182 * t229;
	t246 = t177 * t196;
	t147 = t215 * t246;
	t219 = t174 * t243 - t175 * t182;
	t231 = t175 * t202 * t203;
	t143 = qJD(2) * t231 + t219 * t147 - t174 * t168;
	t255 = t143 * t155 * t156;
	t198 = qJD(3) + qJD(5);
	t239 = qJD(1) * t202;
	t216 = -t186 * t198 + t206 * t239;
	t167 = t213 * qJD(1) + t214 * qJD(2);
	t224 = t198 * t244 + t167;
	t148 = t224 * t193 - t216 * t194;
	t172 = t186 * t193 - t194 * t244;
	t162 = t172 ^ 2;
	t161 = t162 * t164 + 0.1e1;
	t250 = t164 * t172;
	t149 = t216 * t193 + t224 * t194;
	t252 = t149 * t163 * t164;
	t254 = (t148 * t250 - t162 * t252) / t161 ^ 2;
	t245 = t200 * t203;
	t218 = t182 * t245 - t199 * t213;
	t150 = t218 * t246;
	t145 = t219 * t150 + t174 * t213 + t231;
	t253 = t145 * t214;
	t201 = t199 * t200;
	t251 = (t168 * t182 * t200 + t180 * t201 * t238) * t197 / t179 ^ 2;
	t222 = qJD(2) * t256 + qJD(1);
	t237 = qJD(2) * t205;
	t166 = -qJD(1) * t225 - t206 * t237 + t222 * t241;
	t249 = t166 * t156;
	t248 = t174 * t214;
	t247 = t175 * t214;
	t242 = t202 * t206;
	t181 = t214 ^ 2;
	t153 = t181 * t156 + 0.1e1;
	t236 = 0.2e1 * (-t181 * t255 + t214 * t249) / t153 ^ 2;
	t235 = 0.2e1 * t255;
	t234 = 0.2e1 * t254;
	t233 = -0.2e1 * t251;
	t232 = t172 * t252;
	t169 = -qJD(1) * t226 - t204 * t238 + t222 * t240;
	t223 = t198 * t242 - t169;
	t220 = t193 * t163 - t194 * t250;
	t217 = t198 * t213 + t204 * t239;
	t171 = t193 * t242 + t194 * t213;
	t170 = t193 * t213 - t194 * t242;
	t159 = 0.1e1 / t161;
	t151 = 0.1e1 / t153;
	t146 = t257 * t214;
	t142 = (t218 * t233 + (t168 * t245 + t169 * t199 + (-t213 * t245 + (0.2e1 * t201 * t203 ^ 2 + t199) * t182) * qJD(2)) * t177) * t196;
	t140 = -0.2e1 * t254 + 0.2e1 * (t148 * t164 * t159 + (-t159 * t252 - t164 * t254) * t172) * t172;
	t1 = [(-t214 * t199 * t233 + (-t166 * t199 - t214 * t229) * t177) * t196, t142, 0, 0, 0, 0; t182 * t155 * t236 + (-t168 * t155 + (t143 * t182 + t146 * t166) * t156) * t151 - ((t146 * t235 - t257 * t249) * t151 + (t146 * t236 + ((t147 * t177 * t230 + t233) * t248 + (0.2e1 * t230 * t251 - t147 + (-t215 * t196 + t147) * t177) * t247) * t151) * t156) * t214, (-t155 * t186 - t156 * t253) * t236 + (-t235 * t253 + t167 * t155 + (-t186 * t143 + t145 * t166 + (t202 * t237 - t142 * t182 - t150 * t168 + (t150 * t243 + t213) * t147) * t247 + (t147 * t150 * t182 - t169 + (t142 * t205 + (-qJD(2) * t150 - t147) * t203) * t202) * t248) * t156) * t151, 0, 0, 0, 0; (-t163 * t170 + t171 * t250) * t234 + ((t223 * t193 + t217 * t194) * t163 + 0.2e1 * t171 * t232 + (-t170 * t149 - (-t217 * t193 + t223 * t194) * t172 - t171 * t148) * t164) * t159, -t220 * t214 * t234 + (t220 * t166 - ((-t163 * t198 - 0.2e1 * t232) * t194 + (t148 * t194 + (-t172 * t198 + t149) * t193) * t164) * t214) * t159, t140, 0, t140, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:13
	% EndTime: 2019-10-10 12:04:15
	% DurationCPUTime: 2.45s
	% Computational Cost: add. (17252->153), mult. (18084->300), div. (989->12), fcn. (22873->13), ass. (0->136)
	t303 = qJ(3) + pkin(12) + qJ(5);
	t301 = sin(t303);
	t306 = cos(pkin(6));
	t310 = cos(qJ(2));
	t382 = sin(qJ(1));
	t339 = t382 * t310;
	t308 = sin(qJ(2));
	t311 = cos(qJ(1));
	t359 = t311 * t308;
	t322 = -t306 * t359 - t339;
	t302 = cos(t303);
	t305 = sin(pkin(6));
	t362 = t305 * t311;
	t343 = t302 * t362;
	t274 = -t301 * t322 + t343;
	t364 = t305 * t308;
	t345 = t301 * t364;
	t284 = -t306 * t302 + t345;
	t255 = atan2(-t274, t284);
	t250 = sin(t255);
	t251 = cos(t255);
	t241 = -t250 * t274 + t251 * t284;
	t239 = 0.1e1 / t241 ^ 2;
	t340 = t382 * t308;
	t332 = t306 * t340;
	t358 = t311 * t310;
	t292 = -t332 + t358;
	t341 = t305 * t382;
	t279 = t292 * t301 - t302 * t341;
	t269 = t279 ^ 2;
	t235 = t269 * t239 + 0.1e1;
	t321 = -t306 * t339 - t359;
	t271 = t322 * qJD(1) + t321 * qJD(2);
	t304 = qJD(3) + qJD(5);
	t328 = t304 * t341 + t271;
	t338 = qJD(1) * t362;
	t365 = t302 * t304;
	t245 = t292 * t365 + t328 * t301 - t302 * t338;
	t375 = t245 * t239;
	t268 = t274 ^ 2;
	t282 = 0.1e1 / t284 ^ 2;
	t254 = t268 * t282 + 0.1e1;
	t252 = 0.1e1 / t254;
	t336 = t382 * qJD(2);
	t273 = -qJD(1) * t332 - t308 * t336 + (qJD(2) * t306 + qJD(1)) * t358;
	t296 = t301 * t362;
	t337 = t382 * qJD(1);
	t331 = t305 * t337;
	t247 = t273 * t301 - t304 * t296 - t302 * t331 - t322 * t365;
	t356 = qJD(2) * t310;
	t324 = t304 * t306 + t305 * t356;
	t344 = t302 * t364;
	t266 = t324 * t301 + t304 * t344;
	t281 = 0.1e1 / t284;
	t369 = t274 * t282;
	t327 = -t247 * t281 + t266 * t369;
	t229 = t327 * t252;
	t329 = -t250 * t284 - t251 * t274;
	t224 = t329 * t229 - t250 * t247 + t251 * t266;
	t238 = 0.1e1 / t241;
	t240 = t238 * t239;
	t380 = t224 * t240;
	t354 = 0.2e1 * (-t269 * t380 + t279 * t375) / t235 ^ 2;
	t386 = t266 * t282;
	t342 = t306 * t358;
	t289 = -t340 + t342;
	t363 = t305 * t310;
	t323 = -t281 * t289 + t363 * t369;
	t385 = t301 * t323;
	t248 = (t304 * t322 + t331) * t301 + t273 * t302 - t304 * t343;
	t280 = t292 * t302 + t301 * t341;
	t309 = cos(qJ(6));
	t307 = sin(qJ(6));
	t367 = t321 * t307;
	t263 = t280 * t309 - t367;
	t257 = 0.1e1 / t263;
	t258 = 0.1e1 / t263 ^ 2;
	t384 = -0.2e1 * t274;
	t383 = 0.2e1 * t279;
	t246 = t328 * t302 + (-t292 * t304 + t338) * t301;
	t270 = -qJD(1) * t342 - t311 * t356 + (t306 * t336 + t337) * t308;
	t236 = t263 * qJD(6) + t246 * t307 + t270 * t309;
	t366 = t321 * t309;
	t262 = t280 * t307 + t366;
	t256 = t262 ^ 2;
	t244 = t256 * t258 + 0.1e1;
	t372 = t258 * t262;
	t355 = qJD(6) * t262;
	t237 = t246 * t309 - t270 * t307 - t355;
	t377 = t237 * t257 * t258;
	t379 = (t236 * t372 - t256 * t377) / t244 ^ 2;
	t371 = t281 * t386;
	t378 = (t247 * t369 - t268 * t371) / t254 ^ 2;
	t376 = t239 * t279;
	t374 = t250 * t279;
	t373 = t251 * t279;
	t370 = t274 * t281;
	t368 = t321 * t301;
	t361 = t307 * t257;
	t360 = t309 * t262;
	t357 = qJD(2) * t308;
	t353 = -0.2e1 * t379;
	t352 = 0.2e1 * t379;
	t351 = -0.2e1 * t378;
	t350 = t240 * t383;
	t349 = t281 * t378;
	t348 = t239 * t374;
	t347 = t239 * t373;
	t346 = t262 * t377;
	t335 = 0.2e1 * t346;
	t334 = t371 * t384;
	t276 = -t302 * t322 - t296;
	t330 = -qJD(6) * t302 * t321 + t271;
	t261 = -t276 * t309 + t289 * t307;
	t260 = -t276 * t307 - t289 * t309;
	t326 = t258 * t360 - t361;
	t285 = t306 * t301 + t344;
	t325 = -t276 * t281 + t285 * t369;
	t319 = -t250 + (t251 * t370 + t250) * t252;
	t318 = qJD(6) * t292 + t270 * t302 - t304 * t368;
	t272 = t321 * qJD(1) + t322 * qJD(2);
	t267 = t324 * t302 - t304 * t345;
	t265 = t292 * t307 + t302 * t366;
	t264 = -t292 * t309 + t302 * t367;
	t242 = 0.1e1 / t244;
	t233 = 0.1e1 / t235;
	t232 = t252 * t385;
	t230 = t325 * t252;
	t228 = t319 * t279;
	t226 = (-t250 * t289 + t251 * t363) * t301 + t329 * t232;
	t225 = t329 * t230 - t250 * t276 + t251 * t285;
	t222 = t325 * t351 + (t285 * t334 - t248 * t281 + (t247 * t285 + t266 * t276 + t267 * t274) * t282) * t252;
	t221 = t351 * t385 + (t323 * t365 + (t334 * t363 - t272 * t281 + (t266 * t289 + (t247 * t310 - t274 * t357) * t305) * t282) * t301) * t252;
	t220 = t326 * t279 * t353 + (t326 * t245 + ((-qJD(6) * t257 - 0.2e1 * t346) * t309 + (t236 * t309 + (t237 - t355) * t307) * t258) * t279) * t242;
	t219 = (t225 * t376 - t238 * t280) * t354 + (t225 * t224 * t350 + t246 * t238 + (-t280 * t224 - t225 * t245 - (-t222 * t274 - t230 * t247 + t267 + (-t230 * t284 - t276) * t229) * t373 - (-t222 * t284 - t230 * t266 - t248 + (t230 * t274 - t285) * t229) * t374) * t239) * t233;
	t1 = [t349 * t383 + (-t245 * t281 + t279 * t386) * t252, t221, t222, 0, t222, 0; t274 * t238 * t354 + (-t247 * t238 + (t224 * t274 - t228 * t245) * t239) * t233 + (t228 * t239 * t354 + (0.2e1 * t228 * t380 - (-t229 * t252 * t370 + t351) * t348 - (t349 * t384 - t229 + (t229 - t327) * t252) * t347 - t319 * t375) * t233) * t279, (t226 * t376 - t238 * t368) * t354 + (-t226 * t375 + (t270 * t301 + t321 * t365) * t238 + (t226 * t350 - t239 * t368) * t224 - (-t221 * t274 - t232 * t247 + (-t301 * t357 + t310 * t365) * t305 + (-t232 * t284 - t289 * t301) * t229) * t347 - (-t289 * t365 - t221 * t284 - t232 * t266 - t272 * t301 + (t232 * t274 - t301 * t363) * t229) * t348) * t233, t219, 0, t219, 0; (-t257 * t260 + t261 * t372) * t352 + ((t261 * qJD(6) - t248 * t307 - t272 * t309) * t257 + t261 * t335 + (-t260 * t237 - (-t260 * qJD(6) - t248 * t309 + t272 * t307) * t262 - t261 * t236) * t258) * t242, (-t257 * t264 + t265 * t372) * t352 + (t265 * t335 - t330 * t257 * t309 + t318 * t361 + (-t330 * t262 * t307 - t265 * t236 - t264 * t237 - t318 * t360) * t258) * t242, t220, 0, t220, t353 + 0.2e1 * (t236 * t258 * t242 + (-t242 * t377 - t258 * t379) * t262) * t262;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end