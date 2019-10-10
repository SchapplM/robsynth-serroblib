% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR9
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
%   Wie in S6RRPRRR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:03
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR9_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR9_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.42s
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
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:14
	% DurationCPUTime: 0.83s
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
	t138 = sin(pkin(12));
	t140 = cos(pkin(12));
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
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:14
	% DurationCPUTime: 1.09s
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
	t173 = pkin(12) + qJ(4);
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
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:14
	% DurationCPUTime: 0.98s
	% Computational Cost: add. (2443->93), mult. (4665->202), div. (686->14), fcn. (5929->11), ass. (0->95)
	t204 = cos(qJ(2));
	t205 = cos(qJ(1));
	t255 = cos(pkin(6));
	t226 = t205 * t255;
	t224 = t204 * t226;
	t202 = sin(qJ(2));
	t203 = sin(qJ(1));
	t240 = t203 * t202;
	t181 = -t224 + t240;
	t201 = sin(pkin(6));
	t242 = t201 * t204;
	t175 = atan2(-t181, -t242);
	t173 = sin(t175);
	t174 = cos(t175);
	t179 = t181 ^ 2;
	t196 = 0.1e1 / t201 ^ 2;
	t199 = 0.1e1 / t204 ^ 2;
	t178 = t179 * t196 * t199 + 0.1e1;
	t176 = 0.1e1 / t178;
	t195 = 0.1e1 / t201;
	t198 = 0.1e1 / t204;
	t229 = t181 * t195 * t198;
	t256 = t176 * (t174 * t229 - t173) + t173;
	t157 = -t173 * t181 - t174 * t242;
	t154 = 0.1e1 / t157;
	t227 = t203 * t255;
	t225 = t202 * t227;
	t239 = t205 * t204;
	t185 = -t225 + t239;
	t194 = pkin(12) + qJ(4) + qJ(5);
	t192 = sin(t194);
	t193 = cos(t194);
	t243 = t201 * t203;
	t172 = t185 * t193 + t192 * t243;
	t162 = 0.1e1 / t172;
	t155 = 0.1e1 / t157 ^ 2;
	t163 = 0.1e1 / t172 ^ 2;
	t212 = -t202 * t226 - t203 * t204;
	t213 = -t205 * t202 - t204 * t227;
	t167 = -qJD(1) * t213 - qJD(2) * t212;
	t237 = qJD(2) * t202;
	t228 = t199 * t237;
	t214 = t167 * t198 + t181 * t228;
	t245 = t176 * t195;
	t146 = t214 * t245;
	t218 = t173 * t242 - t174 * t181;
	t230 = t174 * t201 * t202;
	t142 = qJD(2) * t230 + t146 * t218 - t173 * t167;
	t254 = t142 * t154 * t155;
	t197 = qJD(4) + qJD(5);
	t238 = qJD(1) * t201;
	t215 = -t185 * t197 + t205 * t238;
	t166 = qJD(1) * t212 + qJD(2) * t213;
	t223 = t197 * t243 + t166;
	t147 = t192 * t223 - t193 * t215;
	t171 = t185 * t192 - t193 * t243;
	t161 = t171 ^ 2;
	t160 = t161 * t163 + 0.1e1;
	t249 = t163 * t171;
	t148 = t192 * t215 + t193 * t223;
	t251 = t148 * t162 * t163;
	t253 = (t147 * t249 - t161 * t251) / t160 ^ 2;
	t244 = t199 * t202;
	t217 = t181 * t244 - t198 * t212;
	t149 = t217 * t245;
	t144 = t149 * t218 + t173 * t212 + t230;
	t252 = t144 * t213;
	t200 = t198 * t199;
	t250 = (t167 * t181 * t199 + t179 * t200 * t237) * t196 / t178 ^ 2;
	t221 = qJD(2) * t255 + qJD(1);
	t236 = qJD(2) * t204;
	t165 = -qJD(1) * t224 - t205 * t236 + t221 * t240;
	t248 = t165 * t155;
	t247 = t173 * t213;
	t246 = t174 * t213;
	t241 = t201 * t205;
	t180 = t213 ^ 2;
	t152 = t180 * t155 + 0.1e1;
	t235 = 0.2e1 * (-t180 * t254 + t213 * t248) / t152 ^ 2;
	t234 = 0.2e1 * t254;
	t233 = 0.2e1 * t253;
	t232 = -0.2e1 * t250;
	t231 = t171 * t251;
	t168 = -qJD(1) * t225 - t203 * t237 + t221 * t239;
	t222 = t197 * t241 - t168;
	t219 = t192 * t162 - t193 * t249;
	t216 = t197 * t212 + t203 * t238;
	t170 = t192 * t241 + t193 * t212;
	t169 = t192 * t212 - t193 * t241;
	t158 = 0.1e1 / t160;
	t150 = 0.1e1 / t152;
	t145 = t256 * t213;
	t141 = (t217 * t232 + (t167 * t244 + t168 * t198 + (-t212 * t244 + (0.2e1 * t200 * t202 ^ 2 + t198) * t181) * qJD(2)) * t176) * t195;
	t139 = -0.2e1 * t253 + 0.2e1 * (t147 * t163 * t158 + (-t158 * t251 - t163 * t253) * t171) * t171;
	t1 = [(-t213 * t198 * t232 + (-t165 * t198 - t213 * t228) * t176) * t195, t141, 0, 0, 0, 0; t181 * t154 * t235 + (-t167 * t154 + (t142 * t181 + t145 * t165) * t155) * t150 - ((t145 * t234 - t256 * t248) * t150 + (t145 * t235 + ((t146 * t176 * t229 + t232) * t247 + (0.2e1 * t229 * t250 - t146 + (-t195 * t214 + t146) * t176) * t246) * t150) * t155) * t213, (-t154 * t185 - t155 * t252) * t235 + (-t234 * t252 + t166 * t154 + (-t185 * t142 + t144 * t165 + (t201 * t236 - t141 * t181 - t149 * t167 + (t149 * t242 + t212) * t146) * t246 + (t146 * t149 * t181 - t168 + (t141 * t204 + (-qJD(2) * t149 - t146) * t202) * t201) * t247) * t155) * t150, 0, 0, 0, 0; (-t162 * t169 + t170 * t249) * t233 + ((t192 * t222 + t193 * t216) * t162 + 0.2e1 * t170 * t231 + (-t169 * t148 - (-t192 * t216 + t193 * t222) * t171 - t170 * t147) * t163) * t158, -t219 * t213 * t233 + (t219 * t165 - ((-t162 * t197 - 0.2e1 * t231) * t193 + (t147 * t193 + (-t171 * t197 + t148) * t192) * t163) * t213) * t158, 0, t139, t139, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:14
	% EndTime: 2019-10-10 11:03:16
	% DurationCPUTime: 2.53s
	% Computational Cost: add. (17252->153), mult. (18084->300), div. (989->12), fcn. (22873->13), ass. (0->136)
	t302 = pkin(12) + qJ(4) + qJ(5);
	t300 = sin(t302);
	t305 = cos(pkin(6));
	t309 = cos(qJ(2));
	t381 = sin(qJ(1));
	t338 = t381 * t309;
	t307 = sin(qJ(2));
	t310 = cos(qJ(1));
	t358 = t310 * t307;
	t321 = -t305 * t358 - t338;
	t301 = cos(t302);
	t304 = sin(pkin(6));
	t361 = t304 * t310;
	t342 = t301 * t361;
	t273 = -t300 * t321 + t342;
	t363 = t304 * t307;
	t344 = t300 * t363;
	t283 = -t305 * t301 + t344;
	t254 = atan2(-t273, t283);
	t249 = sin(t254);
	t250 = cos(t254);
	t240 = -t249 * t273 + t250 * t283;
	t238 = 0.1e1 / t240 ^ 2;
	t339 = t381 * t307;
	t331 = t305 * t339;
	t357 = t310 * t309;
	t291 = -t331 + t357;
	t340 = t304 * t381;
	t278 = t291 * t300 - t301 * t340;
	t268 = t278 ^ 2;
	t234 = t268 * t238 + 0.1e1;
	t320 = -t305 * t338 - t358;
	t270 = t321 * qJD(1) + t320 * qJD(2);
	t303 = qJD(4) + qJD(5);
	t327 = t303 * t340 + t270;
	t337 = qJD(1) * t361;
	t364 = t301 * t303;
	t244 = t291 * t364 + t327 * t300 - t301 * t337;
	t374 = t244 * t238;
	t267 = t273 ^ 2;
	t281 = 0.1e1 / t283 ^ 2;
	t253 = t267 * t281 + 0.1e1;
	t251 = 0.1e1 / t253;
	t335 = t381 * qJD(2);
	t272 = -qJD(1) * t331 - t307 * t335 + (qJD(2) * t305 + qJD(1)) * t357;
	t295 = t300 * t361;
	t336 = t381 * qJD(1);
	t330 = t304 * t336;
	t246 = t272 * t300 - t303 * t295 - t301 * t330 - t321 * t364;
	t355 = qJD(2) * t309;
	t323 = t303 * t305 + t304 * t355;
	t343 = t301 * t363;
	t265 = t323 * t300 + t303 * t343;
	t280 = 0.1e1 / t283;
	t368 = t273 * t281;
	t326 = -t246 * t280 + t265 * t368;
	t228 = t326 * t251;
	t328 = -t249 * t283 - t250 * t273;
	t223 = t328 * t228 - t249 * t246 + t250 * t265;
	t237 = 0.1e1 / t240;
	t239 = t237 * t238;
	t379 = t223 * t239;
	t353 = 0.2e1 * (-t268 * t379 + t278 * t374) / t234 ^ 2;
	t385 = t265 * t281;
	t341 = t305 * t357;
	t288 = -t339 + t341;
	t362 = t304 * t309;
	t322 = -t280 * t288 + t362 * t368;
	t384 = t300 * t322;
	t247 = (t303 * t321 + t330) * t300 + t272 * t301 - t303 * t342;
	t279 = t291 * t301 + t300 * t340;
	t308 = cos(qJ(6));
	t306 = sin(qJ(6));
	t366 = t320 * t306;
	t262 = t279 * t308 - t366;
	t256 = 0.1e1 / t262;
	t257 = 0.1e1 / t262 ^ 2;
	t383 = -0.2e1 * t273;
	t382 = 0.2e1 * t278;
	t245 = t327 * t301 + (-t291 * t303 + t337) * t300;
	t269 = -qJD(1) * t341 - t310 * t355 + (t305 * t335 + t336) * t307;
	t235 = t262 * qJD(6) + t245 * t306 + t269 * t308;
	t365 = t320 * t308;
	t261 = t279 * t306 + t365;
	t255 = t261 ^ 2;
	t243 = t255 * t257 + 0.1e1;
	t371 = t257 * t261;
	t354 = qJD(6) * t261;
	t236 = t245 * t308 - t269 * t306 - t354;
	t376 = t236 * t256 * t257;
	t378 = (t235 * t371 - t255 * t376) / t243 ^ 2;
	t370 = t280 * t385;
	t377 = (t246 * t368 - t267 * t370) / t253 ^ 2;
	t375 = t238 * t278;
	t373 = t249 * t278;
	t372 = t250 * t278;
	t369 = t273 * t280;
	t367 = t320 * t300;
	t360 = t306 * t256;
	t359 = t308 * t261;
	t356 = qJD(2) * t307;
	t352 = -0.2e1 * t378;
	t351 = 0.2e1 * t378;
	t350 = -0.2e1 * t377;
	t349 = t239 * t382;
	t348 = t280 * t377;
	t347 = t238 * t373;
	t346 = t238 * t372;
	t345 = t261 * t376;
	t334 = 0.2e1 * t345;
	t333 = t370 * t383;
	t275 = -t301 * t321 - t295;
	t329 = -qJD(6) * t301 * t320 + t270;
	t260 = -t275 * t308 + t288 * t306;
	t259 = -t275 * t306 - t288 * t308;
	t325 = t257 * t359 - t360;
	t284 = t305 * t300 + t343;
	t324 = -t275 * t280 + t284 * t368;
	t318 = -t249 + (t250 * t369 + t249) * t251;
	t317 = qJD(6) * t291 + t269 * t301 - t303 * t367;
	t271 = t320 * qJD(1) + t321 * qJD(2);
	t266 = t323 * t301 - t303 * t344;
	t264 = t291 * t306 + t301 * t365;
	t263 = -t291 * t308 + t301 * t366;
	t241 = 0.1e1 / t243;
	t232 = 0.1e1 / t234;
	t231 = t251 * t384;
	t229 = t324 * t251;
	t227 = t318 * t278;
	t225 = (-t249 * t288 + t250 * t362) * t300 + t328 * t231;
	t224 = t328 * t229 - t249 * t275 + t250 * t284;
	t221 = t324 * t350 + (t284 * t333 - t247 * t280 + (t246 * t284 + t265 * t275 + t266 * t273) * t281) * t251;
	t220 = t350 * t384 + (t322 * t364 + (t333 * t362 - t271 * t280 + (t265 * t288 + (t246 * t309 - t273 * t356) * t304) * t281) * t300) * t251;
	t219 = t325 * t278 * t352 + (t325 * t244 + ((-qJD(6) * t256 - 0.2e1 * t345) * t308 + (t235 * t308 + (t236 - t354) * t306) * t257) * t278) * t241;
	t218 = (t224 * t375 - t237 * t279) * t353 + (t224 * t223 * t349 + t245 * t237 + (-t279 * t223 - t224 * t244 - (-t221 * t273 - t229 * t246 + t266 + (-t229 * t283 - t275) * t228) * t372 - (-t221 * t283 - t229 * t265 - t247 + (t229 * t273 - t284) * t228) * t373) * t238) * t232;
	t1 = [t348 * t382 + (-t244 * t280 + t278 * t385) * t251, t220, 0, t221, t221, 0; t273 * t237 * t353 + (-t246 * t237 + (t223 * t273 - t227 * t244) * t238) * t232 + (t227 * t238 * t353 + (0.2e1 * t227 * t379 - (-t228 * t251 * t369 + t350) * t347 - (t348 * t383 - t228 + (t228 - t326) * t251) * t346 - t318 * t374) * t232) * t278, (t225 * t375 - t237 * t367) * t353 + (-t225 * t374 + (t269 * t300 + t320 * t364) * t237 + (t225 * t349 - t238 * t367) * t223 - (-t220 * t273 - t231 * t246 + (-t300 * t356 + t309 * t364) * t304 + (-t231 * t283 - t288 * t300) * t228) * t346 - (-t288 * t364 - t220 * t283 - t231 * t265 - t271 * t300 + (t231 * t273 - t300 * t362) * t228) * t347) * t232, 0, t218, t218, 0; (-t256 * t259 + t260 * t371) * t351 + ((t260 * qJD(6) - t247 * t306 - t271 * t308) * t256 + t260 * t334 + (-t259 * t236 - (-t259 * qJD(6) - t247 * t308 + t271 * t306) * t261 - t260 * t235) * t257) * t241, (-t256 * t263 + t264 * t371) * t351 + (t264 * t334 - t329 * t256 * t308 + t317 * t360 + (-t329 * t261 * t306 - t264 * t235 - t263 * t236 - t317 * t359) * t257) * t241, 0, t219, t219, t352 + 0.2e1 * (t235 * t257 * t241 + (-t241 * t376 - t257 * t378) * t261) * t261;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end