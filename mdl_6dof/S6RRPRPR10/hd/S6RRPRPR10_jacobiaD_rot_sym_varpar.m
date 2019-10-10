% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR10
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
%   Wie in S6RRPRPR10_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:20
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR10_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR10_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:53
	% DurationCPUTime: 0.41s
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
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:53
	% DurationCPUTime: 0.82s
	% Computational Cost: add. (1295->85), mult. (3755->195), div. (644->14), fcn. (4884->11), ass. (0->88)
	t141 = sin(qJ(2));
	t143 = cos(qJ(2));
	t142 = sin(qJ(1));
	t181 = cos(pkin(6));
	t159 = t142 * t181;
	t189 = cos(qJ(1));
	t151 = -t189 * t141 - t143 * t159;
	t156 = t181 * t189;
	t155 = t143 * t156;
	t173 = t142 * t141;
	t123 = -t155 + t173;
	t139 = sin(pkin(6));
	t174 = t139 * t143;
	t117 = atan2(-t123, -t174);
	t115 = sin(t117);
	t116 = cos(t117);
	t97 = -t115 * t123 - t116 * t174;
	t95 = 0.1e1 / t97 ^ 2;
	t182 = t151 * t95;
	t121 = t123 ^ 2;
	t134 = 0.1e1 / t139 ^ 2;
	t136 = 0.1e1 / t143 ^ 2;
	t120 = t121 * t134 * t136 + 0.1e1;
	t118 = 0.1e1 / t120;
	t133 = 0.1e1 / t139;
	t135 = 0.1e1 / t143;
	t165 = t123 * t133 * t135;
	t190 = (t116 * t165 - t115) * t118 + t115;
	t94 = 0.1e1 / t97;
	t157 = t141 * t159;
	t163 = t189 * t143;
	t127 = t163 - t157;
	t138 = sin(pkin(11));
	t140 = cos(pkin(11));
	t175 = t139 * t142;
	t114 = t127 * t140 + t138 * t175;
	t108 = 0.1e1 / t114;
	t109 = 0.1e1 / t114 ^ 2;
	t103 = -qJD(1) * t155 - qJD(2) * t163 + (qJD(2) * t181 + qJD(1)) * t173;
	t122 = t151 ^ 2;
	t150 = -t141 * t156 - t142 * t143;
	t105 = -t151 * qJD(1) - t150 * qJD(2);
	t172 = qJD(2) * t141;
	t161 = t136 * t172;
	t152 = t105 * t135 + t123 * t161;
	t178 = t118 * t133;
	t88 = t152 * t178;
	t84 = (-t123 * t88 + t139 * t172) * t116 + (t88 * t174 - t105) * t115;
	t96 = t94 * t95;
	t187 = t84 * t96;
	t92 = t122 * t95 + 0.1e1;
	t188 = (t103 * t182 - t122 * t187) / t92 ^ 2;
	t113 = t127 * t138 - t140 * t175;
	t107 = t113 ^ 2;
	t100 = t107 * t109 + 0.1e1;
	t104 = t150 * qJD(1) + t151 * qJD(2);
	t160 = t189 * qJD(1);
	t158 = t139 * t160;
	t101 = t104 * t138 - t140 * t158;
	t179 = t109 * t113;
	t102 = t104 * t140 + t138 * t158;
	t180 = t102 * t108 * t109;
	t186 = (t101 * t179 - t107 * t180) / t100 ^ 2;
	t185 = t103 * t95;
	t137 = t135 * t136;
	t184 = 0.1e1 / t120 ^ 2 * (t105 * t123 * t136 + t121 * t137 * t172) * t134;
	t177 = t136 * t141;
	t153 = t123 * t177 - t135 * t150;
	t89 = t153 * t178;
	t183 = t123 * t89;
	t176 = t138 * t108;
	t171 = 0.2e1 * t188;
	t170 = 0.2e1 * t186;
	t169 = -0.2e1 * t184;
	t168 = t115 * t182;
	t167 = t116 * t182;
	t166 = t113 * t180;
	t164 = t139 * t189;
	t162 = qJD(1) * t175;
	t112 = t138 * t164 + t140 * t150;
	t111 = t138 * t150 - t140 * t164;
	t106 = -qJD(1) * t157 - t142 * t172 + (qJD(2) * t156 + t160) * t143;
	t98 = 0.1e1 / t100;
	t90 = 0.1e1 / t92;
	t87 = t190 * t151;
	t85 = (t139 * t141 - t183) * t116 + (t89 * t174 + t150) * t115;
	t83 = (t153 * t169 + (t105 * t177 + t106 * t135 + (-t150 * t177 + (0.2e1 * t137 * t141 ^ 2 + t135) * t123) * qJD(2)) * t118) * t133;
	t1 = [(-t151 * t135 * t169 + (-t103 * t135 - t151 * t161) * t118) * t133, t83, 0, 0, 0, 0; -0.2e1 * t188 * t87 * t182 + t123 * t94 * t171 + (-t105 * t94 + (t103 * t87 + t123 * t84) * t95 - (0.2e1 * t187 * t87 + (t118 * t88 * t165 + t169) * t168 + (0.2e1 * t165 * t184 - t88 + (-t152 * t133 + t88) * t118) * t167 - t190 * t185) * t151) * t90, (-t127 * t94 - t85 * t182) * t171 + (t85 * t185 + t104 * t94 + (-0.2e1 * t151 * t85 * t96 - t127 * t95) * t84 + (-t105 * t89 - t123 * t83 + t150 * t88 + (t88 * t89 + qJD(2)) * t174) * t167 + (t88 * t183 - t106 + (t143 * t83 + (-qJD(2) * t89 - t88) * t141) * t139) * t168) * t90, 0, 0, 0, 0; (-t108 * t111 + t112 * t179) * t170 + ((-t106 * t138 + t140 * t162) * t108 + 0.2e1 * t112 * t166 + (-t111 * t102 - (-t106 * t140 - t138 * t162) * t113 - t112 * t101) * t109) * t98, (-t140 * t179 + t176) * t98 * t103 - (-0.2e1 * t140 * t98 * t166 + t170 * t176 + (t102 * t138 * t98 + (t101 * t98 - 0.2e1 * t113 * t186) * t140) * t109) * t151, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:54
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (1788->92), mult. (4303->201), div. (668->14), fcn. (5516->11), ass. (0->92)
	t178 = sin(qJ(2));
	t179 = sin(qJ(1));
	t231 = cos(pkin(6));
	t201 = t179 * t231;
	t199 = t178 * t201;
	t180 = cos(qJ(2));
	t181 = cos(qJ(1));
	t215 = t181 * t180;
	t162 = -t199 + t215;
	t173 = pkin(11) + qJ(4);
	t169 = sin(t173);
	t170 = cos(t173);
	t177 = sin(pkin(6));
	t219 = t177 * t179;
	t191 = -t162 * t169 + t170 * t219;
	t233 = t191 * qJD(4);
	t200 = t181 * t231;
	t198 = t180 * t200;
	t216 = t179 * t178;
	t158 = -t198 + t216;
	t218 = t177 * t180;
	t152 = atan2(-t158, -t218);
	t150 = sin(t152);
	t151 = cos(t152);
	t156 = t158 ^ 2;
	t172 = 0.1e1 / t177 ^ 2;
	t175 = 0.1e1 / t180 ^ 2;
	t155 = t156 * t172 * t175 + 0.1e1;
	t153 = 0.1e1 / t155;
	t171 = 0.1e1 / t177;
	t174 = 0.1e1 / t180;
	t205 = t158 * t171 * t174;
	t232 = (t151 * t205 - t150) * t153 + t150;
	t134 = -t150 * t158 - t151 * t218;
	t131 = 0.1e1 / t134;
	t149 = t162 * t170 + t169 * t219;
	t143 = 0.1e1 / t149;
	t132 = 0.1e1 / t134 ^ 2;
	t144 = 0.1e1 / t149 ^ 2;
	t188 = -t178 * t200 - t179 * t180;
	t189 = -t181 * t178 - t180 * t201;
	t140 = -t189 * qJD(1) - t188 * qJD(2);
	t213 = qJD(2) * t178;
	t202 = t175 * t213;
	t190 = t140 * t174 + t158 * t202;
	t221 = t153 * t171;
	t123 = t190 * t221;
	t194 = t150 * t218 - t151 * t158;
	t206 = t151 * t177 * t178;
	t119 = qJD(2) * t206 + t194 * t123 - t150 * t140;
	t230 = t119 * t131 * t132;
	t220 = t175 * t178;
	t193 = t158 * t220 - t174 * t188;
	t124 = t193 * t221;
	t120 = t194 * t124 + t150 * t188 + t206;
	t229 = t120 * t189;
	t139 = t188 * qJD(1) + t189 * qJD(2);
	t214 = qJD(1) * t177;
	t203 = t181 * t214;
	t128 = t149 * qJD(4) + t139 * t169 - t170 * t203;
	t142 = t191 ^ 2;
	t137 = t142 * t144 + 0.1e1;
	t224 = t144 * t191;
	t129 = t139 * t170 + t169 * t203 + t233;
	t227 = t129 * t143 * t144;
	t228 = (-t128 * t224 - t142 * t227) / t137 ^ 2;
	t176 = t174 * t175;
	t226 = (t140 * t158 * t175 + t156 * t176 * t213) * t172 / t155 ^ 2;
	t197 = qJD(2) * t231 + qJD(1);
	t212 = qJD(2) * t180;
	t138 = -qJD(1) * t198 - t181 * t212 + t197 * t216;
	t225 = t138 * t132;
	t223 = t150 * t189;
	t222 = t151 * t189;
	t217 = t177 * t181;
	t157 = t189 ^ 2;
	t127 = t157 * t132 + 0.1e1;
	t211 = 0.2e1 * (-t157 * t230 + t189 * t225) / t127 ^ 2;
	t210 = 0.2e1 * t230;
	t209 = 0.2e1 * t228;
	t208 = -0.2e1 * t226;
	t207 = t191 * t227;
	t204 = t179 * t214;
	t195 = t169 * t143 + t170 * t224;
	t192 = -t169 * t188 + t170 * t217;
	t147 = t169 * t217 + t170 * t188;
	t141 = -qJD(1) * t199 - t179 * t213 + t197 * t215;
	t135 = 0.1e1 / t137;
	t125 = 0.1e1 / t127;
	t122 = t232 * t189;
	t118 = (t193 * t208 + (t140 * t220 + t141 * t174 + (-t188 * t220 + (0.2e1 * t176 * t178 ^ 2 + t174) * t158) * qJD(2)) * t153) * t171;
	t1 = [(-t189 * t174 * t208 + (-t138 * t174 - t189 * t202) * t153) * t171, t118, 0, 0, 0, 0; t158 * t131 * t211 + (-t140 * t131 + (t119 * t158 + t122 * t138) * t132) * t125 - ((t122 * t210 - t232 * t225) * t125 + (t122 * t211 + ((t123 * t153 * t205 + t208) * t223 + (0.2e1 * t205 * t226 - t123 + (-t190 * t171 + t123) * t153) * t222) * t125) * t132) * t189, (-t131 * t162 - t132 * t229) * t211 + (-t210 * t229 + t139 * t131 + (-t162 * t119 + t120 * t138 + (t177 * t212 - t118 * t158 - t124 * t140 + (t124 * t218 + t188) * t123) * t222 + (t123 * t124 * t158 - t141 + (t118 * t180 + (-qJD(2) * t124 - t123) * t178) * t177) * t223) * t132) * t125, 0, 0, 0, 0; (t143 * t192 - t147 * t224) * t209 + ((t147 * qJD(4) - t141 * t169 + t170 * t204) * t143 - 0.2e1 * t147 * t207 + (t192 * t129 + (t192 * qJD(4) - t141 * t170 - t169 * t204) * t191 - t147 * t128) * t144) * t135, -t195 * t189 * t209 + (t195 * t138 - ((-qJD(4) * t143 + 0.2e1 * t207) * t170 + (t128 * t170 + (t129 + t233) * t169) * t144) * t189) * t135, 0, -0.2e1 * t228 - 0.2e1 * (t128 * t144 * t135 - (-t135 * t227 - t144 * t228) * t191) * t191, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:53
	% EndTime: 2019-10-10 10:20:54
	% DurationCPUTime: 1.73s
	% Computational Cost: add. (7177->123), mult. (11003->252), div. (691->12), fcn. (14028->11), ass. (0->106)
	t207 = cos(pkin(6));
	t209 = sin(qJ(1));
	t210 = cos(qJ(2));
	t266 = cos(qJ(1));
	t231 = t266 * qJD(2);
	t232 = t266 * qJD(1);
	t206 = sin(pkin(6));
	t237 = t206 * t266;
	t208 = sin(qJ(2));
	t250 = t209 * t208;
	t238 = t207 * t250;
	t248 = qJD(2) * t208;
	t271 = -qJD(1) * t238 - t209 * t248 + (t207 * t231 + t232) * t210 - qJD(4) * t237;
	t205 = pkin(11) + qJ(4);
	t203 = sin(t205);
	t204 = cos(t205);
	t236 = t266 * t208;
	t249 = t209 * t210;
	t218 = -t207 * t236 - t249;
	t174 = -t203 * t218 + t204 * t237;
	t253 = t206 * t208;
	t184 = t203 * t253 - t207 * t204;
	t161 = atan2(-t174, t184);
	t156 = sin(t161);
	t157 = cos(t161);
	t151 = -t156 * t174 + t157 * t184;
	t149 = 0.1e1 / t151 ^ 2;
	t235 = t266 * t210;
	t192 = t235 - t238;
	t252 = t206 * t209;
	t222 = -t192 * t203 + t204 * t252;
	t172 = t222 ^ 2;
	t147 = t172 * t149 + 0.1e1;
	t219 = -t207 * t249 - t236;
	t166 = t218 * qJD(1) + t219 * qJD(2);
	t180 = t192 * t204 + t203 * t252;
	t226 = t206 * t232;
	t152 = t180 * qJD(4) + t166 * t203 - t204 * t226;
	t261 = t152 * t149;
	t171 = t174 ^ 2;
	t182 = 0.1e1 / t184 ^ 2;
	t160 = t171 * t182 + 0.1e1;
	t158 = 0.1e1 / t160;
	t234 = qJD(1) * t252;
	t247 = qJD(4) * t204;
	t154 = t271 * t203 - t204 * t234 - t218 * t247;
	t185 = t207 * t203 + t204 * t253;
	t251 = t206 * t210;
	t233 = qJD(2) * t251;
	t169 = t185 * qJD(4) + t203 * t233;
	t181 = 0.1e1 / t184;
	t256 = t174 * t182;
	t224 = -t154 * t181 + t169 * t256;
	t140 = t224 * t158;
	t225 = -t156 * t184 - t157 * t174;
	t136 = t225 * t140 - t156 * t154 + t157 * t169;
	t148 = 0.1e1 / t151;
	t150 = t148 * t149;
	t264 = t136 * t150;
	t246 = 0.2e1 * (-t172 * t264 - t222 * t261) / t147 ^ 2;
	t270 = t169 * t182;
	t228 = t207 * t235;
	t189 = t228 - t250;
	t220 = -t181 * t189 + t251 * t256;
	t269 = t203 * t220;
	t155 = (qJD(4) * t218 + t234) * t203 + t271 * t204;
	t186 = 0.1e1 / t219;
	t187 = 0.1e1 / t219 ^ 2;
	t268 = -0.2e1 * t174;
	t267 = -0.2e1 * t222;
	t258 = t181 * t270;
	t263 = (t154 * t256 - t171 * t258) / t160 ^ 2;
	t262 = t149 * t222;
	t260 = t156 * t222;
	t259 = t157 * t222;
	t257 = t174 * t181;
	t255 = t180 * t187;
	t254 = t219 * t203;
	t245 = -0.2e1 * t263;
	t153 = t222 * qJD(4) + t166 * t204 + t203 * t226;
	t173 = t180 ^ 2;
	t164 = t173 * t187 + 0.1e1;
	t165 = -qJD(1) * t228 - t210 * t231 + (qJD(2) * t207 + qJD(1)) * t250;
	t188 = t186 * t187;
	t244 = 0.2e1 * (-t173 * t188 * t165 + t153 * t255) / t164 ^ 2;
	t243 = t150 * t267;
	t242 = 0.2e1 * t180 * t188;
	t241 = t181 * t263;
	t240 = t149 * t260;
	t239 = t149 * t259;
	t230 = t258 * t268;
	t176 = -t203 * t237 - t204 * t218;
	t223 = -t176 * t181 + t185 * t256;
	t217 = -t156 + (t157 * t257 + t156) * t158;
	t170 = -t184 * qJD(4) + t204 * t233;
	t167 = t219 * qJD(1) + t218 * qJD(2);
	t162 = 0.1e1 / t164;
	t145 = 0.1e1 / t147;
	t144 = t158 * t269;
	t141 = t223 * t158;
	t139 = t217 * t222;
	t138 = (-t156 * t189 + t157 * t251) * t203 + t225 * t144;
	t137 = t225 * t141 - t156 * t176 + t157 * t185;
	t135 = t223 * t245 + (t185 * t230 - t155 * t181 + (t154 * t185 + t169 * t176 + t170 * t174) * t182) * t158;
	t133 = t245 * t269 + (t220 * t247 + (t230 * t251 - t167 * t181 + (t169 * t189 + (t154 * t210 - t174 * t248) * t206) * t182) * t203) * t158;
	t1 = [t241 * t267 + (-t152 * t181 - t222 * t270) * t158, t133, 0, t135, 0, 0; t174 * t148 * t246 + (-t154 * t148 + (t136 * t174 + t139 * t152) * t149) * t145 - (-t139 * t149 * t246 + (-0.2e1 * t139 * t264 + (-t140 * t158 * t257 + t245) * t240 + (t241 * t268 - t140 + (t140 - t224) * t158) * t239 - t217 * t261) * t145) * t222, (-t138 * t262 - t148 * t254) * t246 + (-t138 * t261 + (t165 * t203 + t219 * t247) * t148 + (t138 * t243 - t149 * t254) * t136 + (-t133 * t174 - t144 * t154 + (-t203 * t248 + t210 * t247) * t206 + (-t144 * t184 - t189 * t203) * t140) * t239 + (-t189 * t247 - t133 * t184 - t144 * t169 - t167 * t203 + (t144 * t174 - t203 * t251) * t140) * t240) * t145, 0, (-t137 * t262 - t148 * t180) * t246 + (t137 * t136 * t243 + t153 * t148 + (-t180 * t136 - t137 * t152 + (-t135 * t174 - t141 * t154 + t170 + (-t141 * t184 - t176) * t140) * t259 + (-t135 * t184 - t141 * t169 - t155 + (t141 * t174 - t185) * t140) * t260) * t149) * t145, 0, 0; (-t176 * t186 + t189 * t255) * t244 + (t155 * t186 + t189 * t165 * t242 + (-t189 * t153 - t165 * t176 - t167 * t180) * t187) * t162, (t186 * t204 * t219 + t192 * t255) * t244 + (qJD(4) * t186 * t254 + (-t153 * t192 - t166 * t180) * t187 + (t192 * t242 + (t187 * t219 - t186) * t204) * t165) * t162, 0, t222 * t186 * t244 + (t165 * t187 * t222 + t152 * t186) * t162, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:53
	% EndTime: 2019-10-10 10:20:55
	% DurationCPUTime: 2.09s
	% Computational Cost: add. (8428->151), mult. (13478->298), div. (726->12), fcn. (17045->13), ass. (0->129)
	t270 = sin(qJ(1));
	t273 = cos(qJ(1));
	t267 = cos(pkin(6));
	t290 = qJD(2) * t267 + qJD(1);
	t269 = sin(qJ(2));
	t317 = t270 * t269;
	t299 = t267 * t317;
	t266 = sin(pkin(6));
	t307 = qJD(4) * t266;
	t311 = qJD(2) * t269;
	t272 = cos(qJ(2));
	t313 = t273 * t272;
	t344 = -qJD(1) * t299 - t270 * t311 - t273 * t307 + t290 * t313;
	t253 = -t299 + t313;
	t265 = pkin(11) + qJ(4);
	t263 = sin(t265);
	t264 = cos(t265);
	t321 = t266 * t270;
	t241 = t253 * t263 - t264 * t321;
	t271 = cos(qJ(6));
	t314 = t273 * t269;
	t316 = t270 * t272;
	t252 = t267 * t316 + t314;
	t268 = sin(qJ(6));
	t324 = t252 * t268;
	t286 = t241 * t271 - t324;
	t343 = qJD(6) * t286;
	t251 = t267 * t314 + t316;
	t319 = t266 * t273;
	t237 = t251 * t264 - t263 * t319;
	t322 = t266 * t269;
	t249 = t267 * t263 + t264 * t322;
	t224 = atan2(-t237, t249);
	t211 = sin(t224);
	t212 = cos(t224);
	t202 = -t211 * t237 + t212 * t249;
	t200 = 0.1e1 / t202 ^ 2;
	t242 = t253 * t264 + t263 * t321;
	t235 = t242 ^ 2;
	t196 = t235 * t200 + 0.1e1;
	t229 = -qJD(1) * t251 - qJD(2) * t252;
	t312 = qJD(1) * t266;
	t296 = t273 * t312;
	t309 = qJD(4) * t263;
	t207 = t263 * t296 - t253 * t309 + (t270 * t307 + t229) * t264;
	t332 = t207 * t200;
	t234 = t237 ^ 2;
	t246 = 0.1e1 / t249 ^ 2;
	t223 = t234 * t246 + 0.1e1;
	t217 = 0.1e1 / t223;
	t291 = t344 * t264;
	t297 = t270 * t312;
	t209 = -t251 * t309 + t263 * t297 + t291;
	t248 = -t263 * t322 + t267 * t264;
	t310 = qJD(2) * t272;
	t295 = t266 * t310;
	t233 = qJD(4) * t248 + t264 * t295;
	t245 = 0.1e1 / t249;
	t326 = t237 * t246;
	t285 = -t209 * t245 + t233 * t326;
	t190 = t285 * t217;
	t288 = -t211 * t249 - t212 * t237;
	t185 = t190 * t288 - t211 * t209 + t212 * t233;
	t199 = 0.1e1 / t202;
	t201 = t199 * t200;
	t337 = t185 * t201;
	t306 = 0.2e1 * (-t235 * t337 + t242 * t332) / t196 ^ 2;
	t342 = t233 * t246;
	t298 = t267 * t313;
	t250 = t298 - t317;
	t320 = t266 * t272;
	t282 = -t245 * t250 + t320 * t326;
	t341 = t264 * t282;
	t323 = t252 * t271;
	t222 = t241 * t268 + t323;
	t214 = 0.1e1 / t222;
	t215 = 0.1e1 / t222 ^ 2;
	t340 = -0.2e1 * t237;
	t339 = 0.2e1 * t242;
	t206 = qJD(4) * t242 + t229 * t263 - t264 * t296;
	t228 = -qJD(1) * t298 - t273 * t310 + t290 * t317;
	t197 = qJD(6) * t222 - t206 * t271 - t228 * t268;
	t213 = t286 ^ 2;
	t205 = t213 * t215 + 0.1e1;
	t329 = t215 * t286;
	t198 = t206 * t268 - t228 * t271 + t343;
	t334 = t198 * t214 * t215;
	t336 = (-t197 * t329 - t213 * t334) / t205 ^ 2;
	t328 = t245 * t342;
	t335 = (t209 * t326 - t234 * t328) / t223 ^ 2;
	t333 = t200 * t242;
	t331 = t211 * t242;
	t330 = t212 * t242;
	t327 = t237 * t245;
	t325 = t252 * t264;
	t318 = t268 * t286;
	t315 = t271 * t214;
	t308 = qJD(4) * t264;
	t305 = 0.2e1 * t336;
	t304 = -0.2e1 * t335;
	t303 = t201 * t339;
	t302 = t245 * t335;
	t301 = t200 * t331;
	t300 = t200 * t330;
	t293 = -0.2e1 * t286 * t334;
	t292 = t328 * t340;
	t289 = -qJD(6) * t252 * t263 + t229;
	t236 = t251 * t263 + t264 * t319;
	t287 = -t236 * t271 - t250 * t268;
	t220 = -t236 * t268 + t250 * t271;
	t284 = -t318 * t215 + t315;
	t283 = t236 * t245 + t248 * t326;
	t281 = -t211 + (t212 * t327 + t211) * t217;
	t208 = t251 * t308 + t344 * t263 - t264 * t297;
	t280 = qJD(6) * t253 - t228 * t263 + t252 * t308;
	t232 = -qJD(4) * t249 - t263 * t295;
	t230 = -qJD(1) * t252 - qJD(2) * t251;
	t226 = t253 * t271 - t263 * t324;
	t225 = t253 * t268 + t263 * t323;
	t203 = 0.1e1 / t205;
	t194 = 0.1e1 / t196;
	t193 = t217 * t341;
	t191 = t283 * t217;
	t189 = t281 * t242;
	t187 = (-t211 * t250 + t212 * t320) * t264 + t288 * t193;
	t186 = t191 * t288 + t211 * t236 + t212 * t248;
	t184 = t283 * t304 + (t248 * t292 + t208 * t245 + (t209 * t248 + t232 * t237 - t233 * t236) * t246) * t217;
	t182 = t304 * t341 + (-t282 * t309 + (t292 * t320 - t230 * t245 + (t233 * t250 + (t209 * t272 - t237 * t311) * t266) * t246) * t264) * t217;
	t1 = [t302 * t339 + (-t207 * t245 + t242 * t342) * t217, t182, 0, t184, 0, 0; t237 * t199 * t306 + (((qJD(4) * t251 - t297) * t263 - t291) * t199 + (t185 * t237 - t189 * t207) * t200) * t194 + (t189 * t200 * t306 + (0.2e1 * t189 * t337 - (-t190 * t217 * t327 + t304) * t301 - (t302 * t340 - t190 + (t190 - t285) * t217) * t300 - t281 * t332) * t194) * t242, (t187 * t333 + t199 * t325) * t306 + (-t187 * t332 + (t228 * t264 + t252 * t309) * t199 + (t187 * t303 + t200 * t325) * t185 - (-t182 * t237 - t193 * t209 + (-t264 * t311 - t272 * t309) * t266 + (-t193 * t249 - t250 * t264) * t190) * t300 - (t250 * t309 - t182 * t249 - t193 * t233 - t230 * t264 + (t193 * t237 - t264 * t320) * t190) * t301) * t194, 0, (t186 * t333 + t199 * t241) * t306 + (t186 * t185 * t303 - t206 * t199 + (t241 * t185 - t186 * t207 - (-t184 * t237 - t191 * t209 + t232 + (-t191 * t249 + t236) * t190) * t330 - (-t184 * t249 - t191 * t233 + t208 + (t191 * t237 - t248) * t190) * t331) * t200) * t194, 0, 0; (t214 * t287 - t220 * t329) * t305 + ((qJD(6) * t220 + t208 * t271 + t230 * t268) * t214 + t220 * t293 + (t287 * t198 + (qJD(6) * t287 - t208 * t268 + t230 * t271) * t286 - t220 * t197) * t215) * t203, (-t214 * t225 - t226 * t329) * t305 + (t226 * t293 + t289 * t214 * t268 + t280 * t315 + (t271 * t286 * t289 - t226 * t197 - t225 * t198 - t280 * t318) * t215) * t203, 0, t284 * t242 * t305 + (-t284 * t207 + ((qJD(6) * t214 + t293) * t268 + (-t197 * t268 + (t198 + t343) * t271) * t215) * t242) * t203, 0, -0.2e1 * t336 - 0.2e1 * (t197 * t215 * t203 - (-t203 * t334 - t215 * t336) * t286) * t286;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end