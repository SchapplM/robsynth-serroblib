% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR4
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
%   Wie in S6RRPPRR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:41
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
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
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:14
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (454->46), mult. (1493->122), div. (133->12), fcn. (1882->11), ass. (0->62)
	t133 = sin(qJ(1));
	t134 = cos(qJ(1));
	t131 = cos(pkin(6));
	t129 = sin(pkin(11));
	t132 = sin(qJ(2));
	t161 = cos(pkin(11));
	t166 = cos(qJ(2));
	t145 = -t132 * t129 + t166 * t161;
	t144 = t145 * t131;
	t146 = t129 * t166 + t132 * t161;
	t99 = -t133 * t144 - t134 * t146;
	t92 = t99 ^ 2;
	t109 = t146 * t131;
	t147 = t133 * t109 - t134 * t145;
	t94 = 0.1e1 / t147 ^ 2;
	t169 = t92 * t94;
	t125 = 0.1e1 / t131 ^ 2;
	t130 = sin(pkin(6));
	t123 = t130 ^ 2;
	t128 = t134 ^ 2;
	t119 = t128 * t123 * t125 + 0.1e1;
	t127 = t133 ^ 2;
	t159 = 0.1e1 / t119 ^ 2 * t127;
	t168 = t125 * t159;
	t167 = qJD(1) * t144 + t145 * qJD(2);
	t93 = 0.1e1 / t147;
	t155 = t134 * t130;
	t118 = atan2(t155, t131);
	t114 = sin(t118);
	t115 = cos(t118);
	t105 = t114 * t155 + t115 * t131;
	t102 = 0.1e1 / t105;
	t124 = 0.1e1 / t131;
	t103 = 0.1e1 / t105 ^ 2;
	t163 = t94 * t99;
	t111 = t146 * qJD(2);
	t108 = t131 * t111;
	t156 = t133 * t146;
	t83 = -qJD(1) * t156 - t133 * t108 + t167 * t134;
	t151 = t83 * t163;
	t107 = qJD(2) * t144;
	t97 = -t134 * t109 - t133 * t145;
	t84 = qJD(1) * t97 - t133 * t107 - t134 * t111;
	t95 = t93 * t94;
	t162 = t95 * t84;
	t88 = 0.1e1 + t169;
	t165 = (t162 * t92 - t151) / t88 ^ 2;
	t160 = t103 * t133;
	t158 = t123 * t124;
	t154 = qJD(1) * t134;
	t152 = -0.2e1 * t97 * t99;
	t116 = 0.1e1 / t119;
	t150 = (t116 - 0.1e1) * t130;
	t149 = -0.2e1 * t124 * t168;
	t85 = (-t115 * t116 * t134 * t158 + t114 * t150) * t133;
	t122 = t130 * t123;
	t104 = t102 * t103;
	t96 = t134 * t144 - t156;
	t91 = t127 * t123 * t103 + 0.1e1;
	t86 = 0.1e1 / t88;
	t82 = qJD(1) * t85;
	t1 = [(-t116 * t124 * t130 + t122 * t149) * t154, 0, 0, 0, 0, 0; (0.2e1 * (-t102 * t134 + t160 * t85) / t91 ^ 2 * (-t104 * t127 * t82 + t154 * t160) * t123 + ((0.2e1 * t104 * t133 * t85 - t103 * t134) * t82 + (-t133 * t102 + ((-t85 + (-t122 * t168 - t150) * t133 * t114) * t134 - (t128 * t123 ^ 2 * t149 + (-t159 + (0.2e1 * t127 - t128) * t116) * t158) * t133 * t115) * t103) * qJD(1)) / t91) * t130, 0, 0, 0, 0, 0; (t94 * t152 + 0.2e1 * t96 * t93) * t165 + ((-t97 * t83 - t96 * t84) * t94 - (-t134 * t108 - t167 * t133 - t146 * t154) * t93 + (qJD(1) * t147 - t134 * t107 + t133 * t111) * t163 - t152 * t162) * t86, 0.2e1 * (-t147 * t93 - t169) * t165 + (-0.2e1 * t151 + (t147 * t94 + 0.2e1 * t92 * t95 - t93) * t84) * t86, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:14
	% DurationCPUTime: 1.06s
	% Computational Cost: add. (2643->74), mult. (7979->173), div. (436->14), fcn. (10394->11), ass. (0->83)
	t160 = sin(pkin(6));
	t201 = sin(pkin(11));
	t202 = cos(pkin(11));
	t207 = sin(qJ(2));
	t208 = cos(qJ(2));
	t172 = t208 * t201 + t207 * t202;
	t146 = t172 * t160;
	t137 = qJD(2) * t146;
	t149 = t207 * t201 - t208 * t202;
	t145 = t149 * t160;
	t142 = 0.1e1 / t145 ^ 2;
	t193 = t137 * t142;
	t203 = cos(pkin(6));
	t181 = t203 * t201;
	t182 = t203 * t202;
	t210 = t207 * t181 - t208 * t182;
	t171 = t149 * qJD(2);
	t161 = sin(qJ(1));
	t162 = cos(qJ(1));
	t129 = -t161 * t172 - t162 * t210;
	t121 = atan2(t129, t145);
	t116 = sin(t121);
	t117 = cos(t121);
	t125 = t129 ^ 2;
	t120 = t125 * t142 + 0.1e1;
	t118 = 0.1e1 / t120;
	t141 = 0.1e1 / t145;
	t195 = t129 * t141;
	t209 = (t117 * t195 - t116) * t118 + t116;
	t110 = t116 * t129 + t117 * t145;
	t107 = 0.1e1 / t110;
	t157 = 0.1e1 / t161;
	t108 = 0.1e1 / t110 ^ 2;
	t158 = 0.1e1 / t161 ^ 2;
	t169 = t161 * t210;
	t132 = -t162 * t172 + t169;
	t126 = t132 ^ 2;
	t106 = t126 * t108 + 0.1e1;
	t147 = t208 * t181 + t207 * t182;
	t140 = t147 * qJD(2);
	t112 = t129 * qJD(1) - t161 * t140 - t162 * t171;
	t198 = t112 * t108;
	t190 = qJD(1) * t162;
	t115 = qJD(1) * t169 - t162 * t140 + t161 * t171 - t172 * t190;
	t176 = t115 * t141 - t129 * t193;
	t101 = t176 * t118;
	t180 = -t116 * t145 + t117 * t129;
	t98 = t180 * t101 + t116 * t115 + t117 * t137;
	t205 = t107 * t108 * t98;
	t206 = 0.1e1 / t106 ^ 2 * (-t126 * t205 - t132 * t198);
	t178 = -t162 * t147 + t161 * t149;
	t194 = t129 * t146;
	t175 = -t141 * t178 + t142 * t194;
	t102 = t175 * t118;
	t99 = -t180 * t102 + t116 * t178 + t117 * t146;
	t204 = t132 * t99;
	t192 = t141 * t193;
	t200 = (t129 * t142 * t115 - t125 * t192) / t120 ^ 2;
	t139 = t210 * qJD(2);
	t148 = t172 * qJD(2);
	t113 = t178 * qJD(1) + t161 * t139 - t162 * t148;
	t177 = t161 * t147 + t162 * t149;
	t127 = t177 ^ 2;
	t156 = 0.1e1 / t160 ^ 2;
	t124 = t127 * t158 * t156 + 0.1e1;
	t159 = t157 * t158;
	t199 = (-t113 * t158 * t177 - t127 * t159 * t190) * t156 / t124 ^ 2;
	t197 = t116 * t132;
	t196 = t117 * t132;
	t191 = t158 * t162;
	t189 = -0.2e1 * t206;
	t188 = -0.2e1 * t205;
	t187 = 0.2e1 * t200;
	t186 = -0.2e1 * t141 * t200;
	t179 = t162 * t139 + t161 * t148;
	t155 = 0.1e1 / t160;
	t138 = t160 * t171;
	t122 = 0.1e1 / t124;
	t114 = t177 * qJD(1) + t179;
	t104 = 0.1e1 / t106;
	t100 = t209 * t132;
	t97 = t175 * t187 + (0.2e1 * t192 * t194 + t114 * t141 + (-t115 * t146 + t129 * t138 - t137 * t178) * t142) * t118;
	t1 = [t132 * t186 + (-t112 * t141 - t132 * t193) * t118, t97, 0, 0, 0, 0; t129 * t107 * t189 + (t115 * t107 + (-t100 * t112 - t129 * t98) * t108) * t104 + ((t100 * t188 - t209 * t198) * t104 + (t100 * t189 + ((-t101 * t118 * t195 + t187) * t197 + (t129 * t186 + t101 + (-t101 + t176) * t118) * t196) * t104) * t108) * t132, 0.2e1 * (t107 * t177 - t108 * t204) * t206 + (t113 * t107 + t188 * t204 + (-t99 * t112 + t177 * t98 + (-t102 * t115 + t129 * t97 - t138 + (t102 * t145 + t178) * t101) * t196 + (t102 * t137 - t145 * t97 + t114 + (t102 * t129 - t146) * t101) * t197) * t108) * t104, 0, 0, 0, 0; (0.2e1 * (-t157 * t178 - t177 * t191) * t199 + (t179 * t157 - t113 * t191 + (-0.2e1 * t162 ^ 2 * t159 * t177 - t178 * t191) * qJD(1)) * t122) * t155, (-0.2e1 * t132 * t157 * t199 + (-t132 * t158 * t190 - t112 * t157) * t122) * t155, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:14
	% EndTime: 2019-10-10 09:41:15
	% DurationCPUTime: 1.29s
	% Computational Cost: add. (3306->98), mult. (9869->206), div. (448->12), fcn. (12731->13), ass. (0->96)
	t223 = cos(qJ(1));
	t221 = sin(qJ(1));
	t220 = sin(qJ(2));
	t271 = cos(pkin(11));
	t272 = cos(pkin(6));
	t240 = t272 * t271;
	t217 = sin(pkin(11));
	t244 = t217 * t272;
	t273 = cos(qJ(2));
	t231 = -t220 * t240 - t273 * t244;
	t230 = t221 * t231;
	t233 = -t220 * t217 + t273 * t271;
	t277 = t223 * t233 + t230;
	t211 = -t220 * t244 + t273 * t240;
	t203 = t211 * qJD(2);
	t234 = t273 * t217 + t220 * t271;
	t232 = qJD(2) * t234;
	t276 = t221 * t203 + t223 * t232;
	t218 = sin(pkin(6));
	t209 = t233 * t218;
	t202 = qJD(2) * t209;
	t210 = t234 * t218;
	t206 = 0.1e1 / t210 ^ 2;
	t260 = t202 * t206;
	t193 = t221 * t233 - t223 * t231;
	t180 = atan2(-t193, t210);
	t175 = sin(t180);
	t176 = cos(t180);
	t190 = t193 ^ 2;
	t179 = t190 * t206 + 0.1e1;
	t177 = 0.1e1 / t179;
	t205 = 0.1e1 / t210;
	t262 = t193 * t205;
	t274 = (t176 * t262 + t175) * t177 - t175;
	t164 = -t175 * t193 + t176 * t210;
	t161 = 0.1e1 / t164;
	t196 = -t221 * t211 - t223 * t234;
	t219 = sin(qJ(5));
	t222 = cos(qJ(5));
	t258 = t218 * t221;
	t186 = -t196 * t219 + t222 * t258;
	t182 = 0.1e1 / t186;
	t162 = 0.1e1 / t164 ^ 2;
	t183 = 0.1e1 / t186 ^ 2;
	t253 = qJD(1) * t223;
	t174 = qJD(1) * t230 + t223 * t203 - t221 * t232 + t233 * t253;
	t238 = -t174 * t205 + t193 * t260;
	t155 = t238 * t177;
	t239 = -t175 * t210 - t176 * t193;
	t151 = t239 * t155 - t175 * t174 + t176 * t202;
	t270 = t151 * t161 * t162;
	t192 = t211 * t223 - t221 * t234;
	t204 = t231 * qJD(2);
	t212 = t233 * qJD(2);
	t171 = t192 * qJD(1) + t221 * t204 + t212 * t223;
	t245 = t218 * t253;
	t165 = t186 * qJD(5) - t171 * t222 + t219 * t245;
	t185 = t196 * t222 + t219 * t258;
	t181 = t185 ^ 2;
	t169 = t181 * t183 + 0.1e1;
	t263 = t183 * t185;
	t252 = qJD(5) * t185;
	t166 = t171 * t219 + t222 * t245 - t252;
	t266 = t166 * t182 * t183;
	t269 = (t165 * t263 - t181 * t266) / t169 ^ 2;
	t259 = t205 * t260;
	t268 = (t174 * t193 * t206 - t190 * t259) / t179 ^ 2;
	t267 = t162 * t277;
	t265 = t175 * t277;
	t264 = t176 * t277;
	t261 = t193 * t209;
	t257 = t218 * t223;
	t254 = qJD(1) * t221;
	t191 = t277 ^ 2;
	t160 = t162 * t191 + 0.1e1;
	t172 = t193 * qJD(1) + t276;
	t251 = 0.2e1 * (-t172 * t267 - t191 * t270) / t160 ^ 2;
	t250 = 0.2e1 * t270;
	t249 = 0.2e1 * t269;
	t248 = -0.2e1 * t268;
	t247 = t205 * t268;
	t246 = t218 * t254;
	t243 = 0.2e1 * t185 * t266;
	t237 = t222 * t182 + t219 * t263;
	t236 = -t192 * t205 + t206 * t261;
	t188 = t192 * t219 + t222 * t257;
	t187 = -t192 * t222 + t219 * t257;
	t201 = t218 * t232;
	t173 = t196 * qJD(1) + t204 * t223 - t221 * t212;
	t167 = 0.1e1 / t169;
	t158 = 0.1e1 / t160;
	t156 = t236 * t177;
	t154 = t274 * t277;
	t152 = t239 * t156 - t175 * t192 + t176 * t209;
	t150 = t236 * t248 + (-0.2e1 * t259 * t261 - t173 * t205 + (t174 * t209 + t192 * t202 - t193 * t201) * t206) * t177;
	t1 = [0.2e1 * t277 * t247 + (t172 * t205 + t260 * t277) * t177, t150, 0, 0, 0, 0; t193 * t161 * t251 + (-t174 * t161 + (t151 * t193 + t154 * t172) * t162) * t158 + (t154 * t250 * t158 + (t154 * t251 + (-(-t155 * t177 * t262 + t248) * t265 - (-0.2e1 * t193 * t247 - t155 + (t155 - t238) * t177) * t264 + t274 * t172) * t158) * t162) * t277, (t152 * t267 - t161 * t196) * t251 + (t152 * t277 * t250 - t171 * t161 + (-t196 * t151 + t152 * t172 - (-t150 * t193 - t156 * t174 - t201 + (-t156 * t210 - t192) * t155) * t264 - (-t150 * t210 - t156 * t202 - t173 + (t156 * t193 - t209) * t155) * t265) * t162) * t158, 0, 0, 0, 0; (-t182 * t187 + t188 * t263) * t249 + ((t188 * qJD(5) - t173 * t222 - t219 * t246) * t182 + t188 * t243 + (-t187 * t166 - (-t187 * qJD(5) + t173 * t219 - t222 * t246) * t185 - t188 * t165) * t183) * t167, t237 * t277 * t249 + (-t237 * (t231 * t253 - t233 * t254 - t276) + ((qJD(5) * t182 + t243) * t219 + (-t165 * t219 + (t166 - t252) * t222) * t183) * t277) * t167, 0, 0, -0.2e1 * t269 + 0.2e1 * (t165 * t183 * t167 + (-t167 * t266 - t183 * t269) * t185) * t185, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:14
	% EndTime: 2019-10-10 09:41:17
	% DurationCPUTime: 3.14s
	% Computational Cost: add. (8421->151), mult. (24081->304), div. (726->12), fcn. (31139->15), ass. (0->131)
	t317 = sin(qJ(1));
	t320 = cos(qJ(1));
	t397 = sin(pkin(11));
	t399 = cos(pkin(6));
	t352 = t399 * t397;
	t398 = cos(pkin(11));
	t353 = t399 * t398;
	t400 = sin(qJ(2));
	t401 = cos(qJ(2));
	t334 = t400 * t352 - t401 * t353;
	t338 = t401 * t397 + t400 * t398;
	t287 = -t317 * t338 - t320 * t334;
	t316 = sin(qJ(5));
	t319 = cos(qJ(5));
	t314 = sin(pkin(6));
	t380 = t314 * t320;
	t343 = -t287 * t319 + t316 * t380;
	t276 = t343 ^ 2;
	t308 = t400 * t397 - t401 * t398;
	t302 = t308 * t314;
	t294 = -t302 * t319 + t399 * t316;
	t292 = 0.1e1 / t294 ^ 2;
	t271 = t276 * t292 + 0.1e1;
	t269 = 0.1e1 / t271;
	t336 = qJD(2) * t308;
	t407 = qJD(1) * t334 + t336;
	t304 = t401 * t352 + t400 * t353;
	t408 = -qJD(1) * t338 - qJD(2) * t304;
	t328 = t407 * t317 + t408 * t320;
	t377 = qJD(1) * t314;
	t363 = t317 * t377;
	t364 = t319 * t380;
	t246 = -qJD(5) * t364 + (-t287 * qJD(5) + t363) * t316 + t328 * t319;
	t295 = t302 * t316 + t399 * t319;
	t303 = t338 * t314;
	t297 = qJD(2) * t303;
	t274 = t295 * qJD(5) - t297 * t319;
	t291 = 0.1e1 / t294;
	t383 = t343 * t292;
	t348 = -t246 * t291 - t274 * t383;
	t230 = t348 * t269;
	t272 = atan2(t343, t294);
	t267 = sin(t272);
	t268 = cos(t272);
	t351 = -t267 * t294 + t268 * t343;
	t225 = t351 * t230 - t267 * t246 + t268 * t274;
	t242 = t267 * t343 + t268 * t294;
	t240 = 0.1e1 / t242 ^ 2;
	t411 = t225 * t240;
	t331 = t317 * t334 - t320 * t338;
	t381 = t314 * t317;
	t278 = -t331 * t316 + t319 * t381;
	t315 = sin(qJ(6));
	t318 = cos(qJ(6));
	t349 = t317 * t304 + t308 * t320;
	t255 = t278 * t315 + t318 * t349;
	t410 = 0.2e1 * t255;
	t239 = 0.1e1 / t242;
	t409 = t239 * t411;
	t277 = t316 * t381 + t331 * t319;
	t359 = 0.2e1 * t277 * t409;
	t327 = t408 * t317 - t407 * t320;
	t362 = t320 * t377;
	t249 = t278 * qJD(5) + t316 * t362 - t327 * t319;
	t389 = t249 * t240;
	t406 = -t389 + t359;
	t405 = t274 * t292;
	t286 = t304 * t320 - t317 * t308;
	t345 = t286 * t291 + t303 * t383;
	t404 = t319 * t345;
	t403 = -qJD(6) * t316 * t349 + t327;
	t382 = t349 * t319;
	t402 = -qJD(5) * t382 + t331 * qJD(6);
	t299 = t334 * qJD(2);
	t307 = t338 * qJD(2);
	t266 = t349 * qJD(1) + t299 * t320 + t317 * t307;
	t247 = t343 * qJD(5) - t328 * t316 + t319 * t363;
	t256 = t278 * t318 - t315 * t349;
	t252 = 0.1e1 / t256;
	t253 = 0.1e1 / t256 ^ 2;
	t275 = t277 ^ 2;
	t238 = t240 * t275 + 0.1e1;
	t396 = (-t275 * t409 + t277 * t389) / t238 ^ 2;
	t250 = -t277 * qJD(5) + t327 * t316 + t319 * t362;
	t335 = -t286 * qJD(1) + t317 * t299 - t307 * t320;
	t233 = t256 * qJD(6) + t250 * t315 - t318 * t335;
	t251 = t255 ^ 2;
	t245 = t251 * t253 + 0.1e1;
	t388 = t253 * t255;
	t375 = qJD(6) * t255;
	t234 = t250 * t318 + t315 * t335 - t375;
	t392 = t234 * t252 * t253;
	t395 = (t233 * t388 - t251 * t392) / t245 ^ 2;
	t385 = t291 * t405;
	t393 = (-t246 * t383 - t276 * t385) / t271 ^ 2;
	t391 = t240 * t277;
	t243 = 0.1e1 / t245;
	t390 = t243 * t253;
	t387 = t267 * t277;
	t386 = t268 * t277;
	t384 = t343 * t291;
	t379 = t316 * t315;
	t378 = t316 * t318;
	t376 = qJD(5) * t316;
	t374 = 0.2e1 * t396;
	t373 = -0.2e1 * t395;
	t372 = -0.2e1 * t393;
	t371 = 0.2e1 * t393;
	t369 = t253 * t395;
	t368 = t233 * t390;
	t367 = t255 * t392;
	t366 = t343 * t385;
	t358 = t291 * t371;
	t357 = 0.2e1 * t367;
	t344 = t287 * t316 + t364;
	t258 = -t286 * t315 + t318 * t344;
	t257 = t286 * t318 + t315 * t344;
	t347 = -t315 * t252 + t318 * t388;
	t346 = -t291 * t344 + t295 * t383;
	t339 = -t267 + (-t268 * t384 + t267) * t269;
	t298 = t314 * t336;
	t273 = -t294 * qJD(5) + t297 * t316;
	t262 = t331 * t315 - t349 * t378;
	t236 = 0.1e1 / t238;
	t235 = t269 * t404;
	t232 = t346 * t269;
	t227 = (t267 * t286 - t268 * t303) * t319 + t351 * t235;
	t226 = -t351 * t232 + t267 * t344 + t268 * t295;
	t224 = t346 * t371 + (0.2e1 * t295 * t366 - t247 * t291 + (t246 * t295 - t273 * t343 - t274 * t344) * t292) * t269;
	t222 = t372 * t404 + (-t345 * t376 + (-0.2e1 * t303 * t366 - t266 * t291 + (-t246 * t303 - t274 * t286 - t298 * t343) * t292) * t319) * t269;
	t1 = [t277 * t358 + (-t249 * t291 + t277 * t405) * t269, t222, 0, 0, t224, 0; -0.2e1 * t343 * t239 * t396 + (-t246 * t239 - t343 * t411 - (t339 * t249 + ((t230 * t269 * t384 + t372) * t267 + (t343 * t358 - t230 + (t230 - t348) * t269) * t268) * t277) * t391) * t236 + (t406 * t236 + t391 * t374) * t339 * t277, (t227 * t391 - t239 * t382) * t374 + ((-t319 * t335 - t349 * t376) * t239 + t406 * t227 + (-t382 * t225 - (t303 * t376 + t222 * t343 - t235 * t246 + t298 * t319 + (-t235 * t294 + t286 * t319) * t230) * t386 - (-t286 * t376 - t222 * t294 - t235 * t274 - t266 * t319 + (-t235 * t343 + t303 * t319) * t230) * t387) * t240) * t236, 0, 0, (t226 * t391 - t239 * t278) * t374 + (t226 * t359 + t250 * t239 + (-t278 * t225 - t226 * t249 - (t224 * t343 + t232 * t246 + t273 + (t232 * t294 + t344) * t230) * t386 - (-t224 * t294 + t232 * t274 - t247 + (t232 * t343 - t295) * t230) * t387) * t240) * t236, 0; 0.2e1 * (-t252 * t257 + t258 * t388) * t395 + ((t258 * qJD(6) - t247 * t315 - t266 * t318) * t252 + t258 * t357 + (-t257 * t234 - (-t257 * qJD(6) - t247 * t318 + t266 * t315) * t255 - t258 * t233) * t253) * t243, (t369 * t410 - t368) * t262 + (-t234 * t390 + t252 * t373) * (-t331 * t318 - t349 * t379) + (t262 * t357 + (t379 * t252 - t378 * t388) * t335 + (t403 * t252 - t402 * t388) * t318 + (t402 * t252 + t403 * t388) * t315) * t243, 0, 0, t347 * t277 * t373 + (t347 * t249 + ((-qJD(6) * t252 - 0.2e1 * t367) * t318 + (t233 * t318 + (t234 - t375) * t315) * t253) * t277) * t243, t373 + (t368 + (-t243 * t392 - t369) * t255) * t410;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end