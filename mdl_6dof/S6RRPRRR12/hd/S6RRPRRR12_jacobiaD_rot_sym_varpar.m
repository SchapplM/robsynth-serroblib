% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR12
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
%   Wie in S6RRPRRR12_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:07
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR12_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR12_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:01
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
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
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.80s
	% Computational Cost: add. (1114->72), mult. (3196->173), div. (656->14), fcn. (4222->9), ass. (0->74)
	t125 = cos(qJ(2));
	t126 = cos(qJ(1));
	t156 = cos(pkin(6));
	t141 = t126 * t156;
	t139 = t125 * t141;
	t123 = sin(qJ(2));
	t124 = sin(qJ(1));
	t152 = t124 * t123;
	t104 = -t139 + t152;
	t122 = sin(pkin(6));
	t114 = 0.1e1 / t122;
	t119 = 0.1e1 / t125;
	t144 = t104 * t114 * t119;
	t153 = t122 * t125;
	t94 = atan2(-t104, -t153);
	t92 = sin(t94);
	t93 = cos(t94);
	t101 = t104 ^ 2;
	t115 = 0.1e1 / t122 ^ 2;
	t120 = 0.1e1 / t125 ^ 2;
	t99 = t101 * t115 * t120 + 0.1e1;
	t95 = 0.1e1 / t99;
	t166 = (t93 * t144 - t92) * t95 + t92;
	t87 = -t92 * t104 - t93 * t153;
	t84 = 0.1e1 / t87;
	t116 = 0.1e1 / t124;
	t117 = 0.1e1 / t124 ^ 2;
	t85 = 0.1e1 / t87 ^ 2;
	t149 = qJD(2) * t123;
	t158 = t125 * t92;
	t163 = t104 * t93;
	t143 = t120 * t149;
	t160 = t114 * t95;
	t134 = -t123 * t141 - t124 * t125;
	t142 = t124 * t156;
	t135 = -t126 * t123 - t125 * t142;
	t90 = -t135 * qJD(1) - t134 * qJD(2);
	t77 = (t104 * t143 + t119 * t90) * t160;
	t74 = -t77 * t163 - t92 * t90 + (t93 * t149 + t77 * t158) * t122;
	t165 = t74 * t84 * t85;
	t154 = t120 * t123;
	t136 = t104 * t154 - t119 * t134;
	t78 = t136 * t160;
	t164 = t77 * t78;
	t162 = t135 * t85;
	t161 = t135 * t93;
	t159 = t119 * t95;
	t157 = t92 * t135;
	t155 = t117 * t126;
	t151 = t126 * t125;
	t150 = qJD(1) * t126;
	t102 = t135 ^ 2;
	t81 = t102 * t85 + 0.1e1;
	t138 = qJD(2) * t156 + qJD(1);
	t88 = -qJD(1) * t139 - qJD(2) * t151 + t138 * t152;
	t148 = 0.2e1 * (-t102 * t165 + t88 * t162) / t81 ^ 2;
	t147 = 0.2e1 * t165;
	t121 = t119 * t120;
	t146 = -0.2e1 * (t101 * t121 * t149 + t104 * t120 * t90) * t115 / t99 ^ 2;
	t140 = t123 * t142;
	t108 = -t140 + t151;
	t103 = t108 ^ 2;
	t100 = t103 * t117 * t115 + 0.1e1;
	t118 = t116 * t117;
	t89 = t134 * qJD(1) + t135 * qJD(2);
	t145 = 0.2e1 * (-t103 * t118 * t150 + t108 * t117 * t89) * t115 / t100 ^ 2;
	t133 = t119 * t146 + t95 * t143;
	t97 = 0.1e1 / t100;
	t91 = -qJD(1) * t140 - t124 * t149 + t138 * t151;
	t79 = 0.1e1 / t81;
	t76 = t166 * t135;
	t75 = -t78 * t163 + t92 * t134 + (t123 * t93 + t78 * t158) * t122;
	t73 = (t136 * t146 + (t90 * t154 + t119 * t91 + (-t134 * t154 + (0.2e1 * t121 * t123 ^ 2 + t119) * t104) * qJD(2)) * t95) * t114;
	t1 = [(-t133 * t135 - t88 * t159) * t114, t73, 0, 0, 0, 0; t104 * t84 * t148 + (-t90 * t84 + (t104 * t74 + t76 * t88) * t85) * t79 - (t76 * t147 * t79 + (t76 * t148 + ((t77 * t95 * t144 + t146) * t157 + ((t95 - 0.1e1) * t77 + (-t133 * t104 - t90 * t159) * t114) * t161 - t166 * t88) * t79) * t85) * t135, (-t108 * t84 - t75 * t162) * t148 + (-t75 * t135 * t147 + t89 * t84 + (-t108 * t74 + t75 * t88 + (t104 * t164 - t91) * t157 + (-t104 * t73 + t134 * t77 - t78 * t90) * t161) * t85 + ((-qJD(2) * t78 - t77) * t92 * t123 + (t73 * t92 + (qJD(2) + t164) * t93) * t125) * t122 * t162) * t79, 0, 0, 0, 0; ((t108 * t155 - t116 * t134) * t145 + (-t89 * t155 - t116 * t91 + (-t134 * t155 + (0.2e1 * t118 * t126 ^ 2 + t116) * t108) * qJD(1)) * t97) * t114, (t116 * t88 * t97 - (t117 * t97 * t150 + t116 * t145) * t135) * t114, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:03
	% DurationCPUTime: 0.88s
	% Computational Cost: add. (1334->90), mult. (4303->199), div. (668->14), fcn. (5516->11), ass. (0->91)
	t174 = sin(qJ(2));
	t175 = sin(qJ(1));
	t177 = cos(qJ(2));
	t178 = cos(qJ(1));
	t225 = cos(pkin(6));
	t194 = t178 * t225;
	t157 = t174 * t194 + t175 * t177;
	t172 = sin(pkin(6));
	t213 = t172 * t174;
	t151 = atan2(-t157, t213);
	t147 = sin(t151);
	t148 = cos(t151);
	t154 = t157 ^ 2;
	t168 = 0.1e1 / t172 ^ 2;
	t170 = 0.1e1 / t174 ^ 2;
	t152 = t154 * t168 * t170 + 0.1e1;
	t149 = 0.1e1 / t152;
	t167 = 0.1e1 / t172;
	t169 = 0.1e1 / t174;
	t199 = t157 * t167 * t169;
	t226 = t149 * (t148 * t199 + t147) - t147;
	t131 = -t147 * t157 + t148 * t213;
	t128 = 0.1e1 / t131;
	t195 = t175 * t225;
	t159 = t178 * t174 + t177 * t195;
	t173 = sin(qJ(4));
	t176 = cos(qJ(4));
	t212 = t172 * t175;
	t144 = t159 * t173 + t176 * t212;
	t140 = 0.1e1 / t144;
	t129 = 0.1e1 / t131 ^ 2;
	t141 = 0.1e1 / t144 ^ 2;
	t190 = qJD(2) * t225 + qJD(1);
	t192 = t174 * t195;
	t207 = qJD(2) * t174;
	t209 = t178 * t177;
	t138 = -qJD(1) * t192 - t175 * t207 + t190 * t209;
	t206 = qJD(2) * t177;
	t196 = t170 * t206;
	t185 = -t138 * t169 + t157 * t196;
	t215 = t149 * t167;
	t120 = t185 * t215;
	t187 = -t147 * t213 - t148 * t157;
	t200 = t148 * t172 * t177;
	t116 = qJD(2) * t200 + t120 * t187 - t147 * t138;
	t224 = t116 * t128 * t129;
	t191 = t177 * t194;
	t210 = t175 * t174;
	t156 = -t191 + t210;
	t214 = t170 * t177;
	t186 = t156 * t169 + t157 * t214;
	t121 = t186 * t215;
	t117 = t121 * t187 + t147 * t156 + t200;
	t160 = -t192 + t209;
	t223 = t117 * t160;
	t135 = -qJD(1) * t191 - t178 * t206 + t190 * t210;
	t208 = qJD(1) * t172;
	t197 = t178 * t208;
	t126 = qJD(4) * t144 + t135 * t176 + t173 * t197;
	t143 = -t159 * t176 + t173 * t212;
	t139 = t143 ^ 2;
	t134 = t139 * t141 + 0.1e1;
	t218 = t141 * t143;
	t205 = qJD(4) * t143;
	t127 = -t135 * t173 + t176 * t197 - t205;
	t220 = t127 * t140 * t141;
	t222 = (t126 * t218 - t139 * t220) / t134 ^ 2;
	t171 = t169 * t170;
	t221 = (t138 * t157 * t170 - t154 * t171 * t206) * t168 / t152 ^ 2;
	t136 = qJD(1) * t157 + qJD(2) * t159;
	t219 = t136 * t129;
	t217 = t147 * t160;
	t216 = t148 * t160;
	t211 = t172 * t178;
	t155 = t160 ^ 2;
	t124 = t155 * t129 + 0.1e1;
	t204 = 0.2e1 * (-t155 * t224 - t160 * t219) / t124 ^ 2;
	t203 = 0.2e1 * t224;
	t202 = 0.2e1 * t222;
	t201 = -0.2e1 * t221;
	t198 = t175 * t208;
	t193 = 0.2e1 * t143 * t220;
	t188 = t176 * t140 + t173 * t218;
	t146 = -t156 * t173 + t176 * t211;
	t145 = t156 * t176 + t173 * t211;
	t137 = qJD(1) * t159 + qJD(2) * t157;
	t132 = 0.1e1 / t134;
	t122 = 0.1e1 / t124;
	t119 = t226 * t160;
	t115 = (t186 * t201 + (t138 * t214 + t137 * t169 + (-t156 * t214 + (-0.2e1 * t171 * t177 ^ 2 - t169) * t157) * qJD(2)) * t149) * t167;
	t1 = [(0.2e1 * t160 * t169 * t221 + (t136 * t169 + t160 * t196) * t149) * t167, t115, 0, 0, 0, 0; t157 * t128 * t204 + (-t138 * t128 + (t116 * t157 + t119 * t136) * t129) * t122 + ((t119 * t203 + t226 * t219) * t122 + (t119 * t204 + (-(-t120 * t149 * t199 + t201) * t217 - (t199 * t201 - t120 + (-t167 * t185 + t120) * t149) * t216) * t122) * t129) * t160, (t128 * t159 + t129 * t223) * t204 + (t203 * t223 + t135 * t128 + (t159 * t116 + t117 * t136 - (-t172 * t207 - t115 * t157 - t121 * t138 + (-t121 * t213 + t156) * t120) * t216 - (t120 * t121 * t157 + t137 + (-t115 * t174 + (-qJD(2) * t121 - t120) * t177) * t172) * t217) * t129) * t122, 0, 0, 0, 0; (-t140 * t145 + t146 * t218) * t202 + ((qJD(4) * t146 + t137 * t176 - t173 * t198) * t140 + t146 * t193 + (-t145 * t127 - (-qJD(4) * t145 - t137 * t173 - t176 * t198) * t143 - t146 * t126) * t141) * t132, t188 * t160 * t202 + (t188 * t136 + ((qJD(4) * t140 + t193) * t173 + (-t126 * t173 + (t127 - t205) * t176) * t141) * t160) * t132, 0, -0.2e1 * t222 + 0.2e1 * (t126 * t141 * t132 + (-t132 * t220 - t141 * t222) * t143) * t143, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:03
	% DurationCPUTime: 1.06s
	% Computational Cost: add. (1905->92), mult. (4665->200), div. (686->14), fcn. (5929->11), ass. (0->95)
	t204 = sin(qJ(2));
	t205 = sin(qJ(1));
	t206 = cos(qJ(2));
	t207 = cos(qJ(1));
	t255 = cos(pkin(6));
	t227 = t207 * t255;
	t184 = t204 * t227 + t205 * t206;
	t203 = sin(pkin(6));
	t243 = t203 * t204;
	t178 = atan2(-t184, t243);
	t174 = sin(t178);
	t175 = cos(t178);
	t181 = t184 ^ 2;
	t197 = 0.1e1 / t203 ^ 2;
	t200 = 0.1e1 / t204 ^ 2;
	t179 = t181 * t197 * t200 + 0.1e1;
	t176 = 0.1e1 / t179;
	t196 = 0.1e1 / t203;
	t199 = 0.1e1 / t204;
	t230 = t184 * t196 * t199;
	t256 = (t175 * t230 + t174) * t176 - t174;
	t158 = -t174 * t184 + t175 * t243;
	t155 = 0.1e1 / t158;
	t228 = t205 * t255;
	t186 = t207 * t204 + t206 * t228;
	t202 = qJ(4) + qJ(5);
	t194 = sin(t202);
	t195 = cos(t202);
	t242 = t203 * t205;
	t171 = t186 * t194 + t195 * t242;
	t167 = 0.1e1 / t171;
	t156 = 0.1e1 / t158 ^ 2;
	t168 = 0.1e1 / t171 ^ 2;
	t221 = qJD(2) * t255 + qJD(1);
	t225 = t204 * t228;
	t237 = qJD(2) * t204;
	t239 = t207 * t206;
	t165 = -qJD(1) * t225 - t205 * t237 + t221 * t239;
	t236 = qJD(2) * t206;
	t229 = t200 * t236;
	t214 = -t165 * t199 + t184 * t229;
	t245 = t176 * t196;
	t147 = t214 * t245;
	t218 = -t174 * t243 - t175 * t184;
	t231 = t175 * t203 * t206;
	t143 = qJD(2) * t231 + t218 * t147 - t174 * t165;
	t254 = t143 * t155 * t156;
	t224 = t206 * t227;
	t240 = t205 * t204;
	t183 = -t224 + t240;
	t244 = t200 * t206;
	t217 = t183 * t199 + t184 * t244;
	t148 = t217 * t245;
	t144 = t218 * t148 + t174 * t183 + t231;
	t187 = -t225 + t239;
	t253 = t144 * t187;
	t198 = qJD(4) + qJD(5);
	t238 = qJD(1) * t203;
	t215 = t186 * t198 + t207 * t238;
	t162 = -qJD(1) * t224 - t207 * t236 + t221 * t240;
	t223 = t198 * t242 + t162;
	t149 = t215 * t194 + t223 * t195;
	t170 = -t186 * t195 + t194 * t242;
	t166 = t170 ^ 2;
	t161 = t166 * t168 + 0.1e1;
	t248 = t168 * t170;
	t150 = -t223 * t194 + t215 * t195;
	t251 = t150 * t167 * t168;
	t252 = (t149 * t248 - t166 * t251) / t161 ^ 2;
	t201 = t199 * t200;
	t250 = (t165 * t184 * t200 - t181 * t201 * t236) * t197 / t179 ^ 2;
	t163 = t184 * qJD(1) + t186 * qJD(2);
	t249 = t163 * t156;
	t247 = t174 * t187;
	t246 = t175 * t187;
	t241 = t203 * t207;
	t182 = t187 ^ 2;
	t153 = t182 * t156 + 0.1e1;
	t235 = 0.2e1 * (-t182 * t254 - t187 * t249) / t153 ^ 2;
	t234 = 0.2e1 * t254;
	t233 = 0.2e1 * t252;
	t232 = -0.2e1 * t250;
	t226 = 0.2e1 * t170 * t251;
	t164 = t186 * qJD(1) + qJD(2) * t184;
	t222 = t198 * t241 + t164;
	t219 = t195 * t167 + t194 * t248;
	t216 = -t183 * t198 - t205 * t238;
	t173 = -t183 * t194 + t195 * t241;
	t172 = t183 * t195 + t194 * t241;
	t159 = 0.1e1 / t161;
	t151 = 0.1e1 / t153;
	t146 = t256 * t187;
	t142 = (t217 * t232 + (t165 * t244 + t164 * t199 + (-t183 * t244 + (-0.2e1 * t201 * t206 ^ 2 - t199) * t184) * qJD(2)) * t176) * t196;
	t140 = -0.2e1 * t252 + 0.2e1 * (t149 * t168 * t159 + (-t159 * t251 - t168 * t252) * t170) * t170;
	t1 = [(0.2e1 * t187 * t199 * t250 + (t163 * t199 + t187 * t229) * t176) * t196, t142, 0, 0, 0, 0; t184 * t155 * t235 + (-t165 * t155 + (t143 * t184 + t146 * t163) * t156) * t151 + ((t146 * t234 + t256 * t249) * t151 + (t146 * t235 + (-(-t147 * t176 * t230 + t232) * t247 - (t230 * t232 - t147 + (-t214 * t196 + t147) * t176) * t246) * t151) * t156) * t187, (t155 * t186 + t156 * t253) * t235 + (t234 * t253 + t162 * t155 + (t186 * t143 + t144 * t163 - (-t203 * t237 - t142 * t184 - t148 * t165 + (-t148 * t243 + t183) * t147) * t246 - (t147 * t148 * t184 + t164 + (-t142 * t204 + (-qJD(2) * t148 - t147) * t206) * t203) * t247) * t156) * t151, 0, 0, 0, 0; (-t167 * t172 + t173 * t248) * t233 + ((t216 * t194 + t222 * t195) * t167 + t173 * t226 + (-t172 * t150 - (-t222 * t194 + t216 * t195) * t170 - t173 * t149) * t168) * t159, t219 * t187 * t233 + (t219 * t163 + ((t167 * t198 + t226) * t194 + (-t149 * t194 + (-t170 * t198 + t150) * t195) * t168) * t187) * t159, 0, t140, t140, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:05
	% DurationCPUTime: 3.12s
	% Computational Cost: add. (11947->152), mult. (18084->300), div. (989->12), fcn. (22873->13), ass. (0->132)
	t295 = qJD(4) + qJD(5);
	t299 = sin(qJ(2));
	t300 = sin(qJ(1));
	t302 = cos(qJ(2));
	t303 = cos(qJ(1));
	t382 = cos(pkin(6));
	t333 = t303 * t382;
	t319 = -t299 * t333 - t300 * t302;
	t334 = t300 * t382;
	t320 = t303 * t299 + t302 * t334;
	t297 = sin(pkin(6));
	t358 = t297 * t303;
	t394 = t320 * qJD(1) - t319 * qJD(2) + t295 * t358;
	t284 = t299 * t300 - t302 * t333;
	t296 = qJ(4) + qJ(5);
	t293 = sin(t296);
	t294 = cos(t296);
	t322 = t284 * t294 + t293 * t358;
	t269 = t322 ^ 2;
	t359 = t297 * t302;
	t339 = t294 * t359;
	t281 = t382 * t293 + t339;
	t279 = 0.1e1 / t281 ^ 2;
	t258 = t269 * t279 + 0.1e1;
	t256 = 0.1e1 / t258;
	t355 = qJD(1) * t300;
	t337 = t297 * t355;
	t392 = t394 * t294;
	t241 = (t284 * t295 + t337) * t293 - t392;
	t354 = qJD(2) * t299;
	t317 = -t382 * t295 + t297 * t354;
	t340 = t293 * t359;
	t265 = -t317 * t294 - t295 * t340;
	t278 = 0.1e1 / t281;
	t368 = t322 * t279;
	t326 = -t241 * t278 - t265 * t368;
	t225 = t326 * t256;
	t259 = atan2(t322, t281);
	t250 = sin(t259);
	t251 = cos(t259);
	t328 = -t250 * t281 + t251 * t322;
	t220 = t328 * t225 - t241 * t250 + t251 * t265;
	t237 = t250 * t322 + t251 * t281;
	t235 = 0.1e1 / t237 ^ 2;
	t393 = t220 * t235;
	t362 = t294 * t295;
	t242 = t284 * t362 + t293 * t394 + t294 * t337;
	t360 = t297 * t300;
	t271 = t320 * t293 + t294 * t360;
	t329 = t299 * t334;
	t357 = t303 * t302;
	t286 = -t329 + t357;
	t298 = sin(qJ(6));
	t301 = cos(qJ(6));
	t252 = t271 * t298 - t286 * t301;
	t391 = 0.2e1 * t252;
	t234 = 0.1e1 / t237;
	t390 = t234 * t393;
	t270 = t293 * t360 - t320 * t294;
	t330 = 0.2e1 * t270 * t390;
	t314 = (t382 * qJD(1) + qJD(2)) * t357 - qJD(2) * t329 - t299 * t355;
	t383 = t295 * t360 - t314;
	t384 = qJD(1) * t358 + t320 * t295;
	t244 = t384 * t293 + t383 * t294;
	t374 = t244 * t235;
	t389 = -t374 + t330;
	t388 = t265 * t279;
	t361 = t297 * t299;
	t343 = t322 * t361;
	t321 = -t278 * t319 + t279 * t343;
	t387 = t294 * t321;
	t386 = -t286 * t293 * qJD(6) - t314;
	t385 = -t320 * qJD(6) + t286 * t362;
	t253 = t271 * t301 + t286 * t298;
	t247 = 0.1e1 / t253;
	t248 = 0.1e1 / t253 ^ 2;
	t268 = t270 ^ 2;
	t231 = t235 * t268 + 0.1e1;
	t381 = (-t268 * t390 + t270 * t374) / t231 ^ 2;
	t245 = -t383 * t293 + t384 * t294;
	t266 = t319 * qJD(1) - t320 * qJD(2);
	t232 = t253 * qJD(6) + t245 * t298 - t266 * t301;
	t246 = t252 ^ 2;
	t240 = t246 * t248 + 0.1e1;
	t373 = t248 * t252;
	t352 = qJD(6) * t252;
	t233 = t245 * t301 + t266 * t298 - t352;
	t377 = t233 * t247 * t248;
	t380 = (t232 * t373 - t246 * t377) / t240 ^ 2;
	t370 = t278 * t388;
	t378 = (-t241 * t368 - t269 * t370) / t258 ^ 2;
	t376 = t235 * t270;
	t238 = 0.1e1 / t240;
	t375 = t238 * t248;
	t372 = t250 * t270;
	t371 = t251 * t270;
	t369 = t322 * t278;
	t282 = t382 * t294 - t340;
	t367 = t322 * t282;
	t366 = t286 * t294;
	t365 = t293 * t295;
	t364 = t293 * t298;
	t363 = t293 * t301;
	t353 = qJD(2) * t302;
	t351 = 0.2e1 * t381;
	t350 = -0.2e1 * t380;
	t349 = -0.2e1 * t378;
	t348 = 0.2e1 * t378;
	t346 = t248 * t380;
	t345 = t232 * t375;
	t344 = t252 * t377;
	t332 = t278 * t348;
	t331 = 0.2e1 * t344;
	t323 = -t284 * t293 + t294 * t358;
	t255 = t298 * t319 + t301 * t323;
	t254 = t298 * t323 - t301 * t319;
	t325 = -t298 * t247 + t301 * t373;
	t324 = -t278 * t323 + t279 * t367;
	t318 = -t250 + (-t251 * t369 + t250) * t256;
	t267 = -qJD(1) * t329 - t300 * t354 + (qJD(2) * t382 + qJD(1)) * t357;
	t264 = t317 * t293 - t295 * t339;
	t261 = t286 * t363 - t320 * t298;
	t229 = 0.1e1 / t231;
	t228 = t256 * t387;
	t227 = t324 * t256;
	t222 = (-t250 * t319 - t251 * t361) * t294 + t328 * t228;
	t221 = -t328 * t227 + t250 * t323 + t251 * t282;
	t218 = t324 * t348 + (0.2e1 * t367 * t370 - t242 * t278 + (t241 * t282 - t264 * t322 - t265 * t323) * t279) * t256;
	t217 = t349 * t387 + (-t321 * t365 + (-0.2e1 * t343 * t370 + t267 * t278 + (t265 * t319 + (-t241 * t299 + t322 * t353) * t297) * t279) * t294) * t256;
	t216 = t325 * t270 * t350 + (t325 * t244 + ((-qJD(6) * t247 - 0.2e1 * t344) * t301 + (t232 * t301 + (t233 - t352) * t298) * t248) * t270) * t238;
	t215 = (t221 * t376 - t234 * t271) * t351 + (t221 * t330 + t245 * t234 + (-t271 * t220 - t221 * t244 - (t218 * t322 + t227 * t241 + t264 + (t227 * t281 + t323) * t225) * t371 - (-t218 * t281 + t227 * t265 - t242 + (t227 * t322 - t282) * t225) * t372) * t235) * t229;
	t1 = [t270 * t332 + (-t244 * t278 + t270 * t388) * t256, t217, 0, t218, t218, 0; -0.2e1 * t322 * t234 * t381 + ((-t284 * t365 - t293 * t337 + t392) * t234 - t322 * t393 - (t318 * t244 + ((t225 * t256 * t369 + t349) * t250 + (t322 * t332 - t225 + (t225 - t326) * t256) * t251) * t270) * t376) * t229 + (t389 * t229 + t376 * t351) * t318 * t270, (t222 * t376 + t234 * t366) * t351 + ((-t266 * t294 + t286 * t365) * t234 + t389 * t222 + (t366 * t220 - (t217 * t322 - t228 * t241 + (-t294 * t353 + t299 * t365) * t297 + (-t228 * t281 - t294 * t319) * t225) * t371 - (t319 * t365 - t217 * t281 - t228 * t265 + t267 * t294 + (-t228 * t322 + t294 * t361) * t225) * t372) * t235) * t229, 0, t215, t215, 0; 0.2e1 * (-t247 * t254 + t255 * t373) * t380 + ((t255 * qJD(6) - t242 * t298 + t267 * t301) * t247 + t255 * t331 + (-t254 * t233 - (-t254 * qJD(6) - t242 * t301 - t267 * t298) * t252 - t255 * t232) * t248) * t238, (t346 * t391 - t345) * t261 + (-t233 * t375 + t247 * t350) * (t286 * t364 + t320 * t301) + (t261 * t331 + (t364 * t247 - t363 * t373) * t266 + (-t386 * t247 - t385 * t373) * t301 + (t385 * t247 - t386 * t373) * t298) * t238, 0, t216, t216, t350 + (t345 + (-t238 * t377 - t346) * t252) * t391;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end