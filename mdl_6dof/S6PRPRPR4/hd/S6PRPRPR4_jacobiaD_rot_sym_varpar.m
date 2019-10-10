% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR4
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
%   Wie in S6PRPRPR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(10));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(10)) * cos(pkin(6));
	t37 = t39 * t42 - t41 * t45;
	t34 = 0.1e1 / t37 ^ 2;
	t49 = qJD(2) * t34;
	t36 = t39 * t41 + t42 * t45;
	t33 = t36 ^ 2;
	t30 = t33 * t34 + 0.1e1;
	t46 = t37 * t49;
	t47 = t36 / t37 * t49;
	t48 = (t33 * t47 + t36 * t46) / t30 ^ 2;
	t28 = 0.1e1 / t30;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, -0.2e1 * t48 + 0.2e1 * (t28 * t46 + (t28 * t47 - t34 * t48) * t36) * t36, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:32
	% DurationCPUTime: 0.52s
	% Computational Cost: add. (610->48), mult. (1820->126), div. (396->14), fcn. (2421->11), ass. (0->64)
	t116 = sin(qJ(2));
	t117 = cos(qJ(2));
	t113 = sin(pkin(10));
	t136 = cos(pkin(6));
	t129 = t113 * t136;
	t135 = cos(pkin(10));
	t104 = -t116 * t129 + t135 * t117;
	t126 = t136 * t135;
	t100 = t113 * t116 - t117 * t126;
	t114 = sin(pkin(6));
	t132 = t114 * t117;
	t90 = atan2(-t100, -t132);
	t88 = sin(t90);
	t89 = cos(t90);
	t77 = -t88 * t100 - t89 * t132;
	t74 = 0.1e1 / t77;
	t112 = sin(pkin(11));
	t115 = cos(pkin(11));
	t133 = t113 * t114;
	t87 = t104 * t115 + t112 * t133;
	t83 = 0.1e1 / t87;
	t109 = 0.1e1 / t117;
	t110 = 0.1e1 / t117 ^ 2;
	t75 = 0.1e1 / t77 ^ 2;
	t84 = 0.1e1 / t87 ^ 2;
	t131 = qJD(2) * t116;
	t137 = t117 * t88;
	t143 = t100 * t89;
	t134 = t110 * t116;
	t130 = t100 * t134;
	t107 = 0.1e1 / t114;
	t108 = 0.1e1 / t114 ^ 2;
	t98 = t100 ^ 2;
	t93 = t98 * t108 * t110 + 0.1e1;
	t91 = 0.1e1 / t93;
	t141 = t107 * t91;
	t102 = t113 * t117 + t116 * t126;
	t95 = t102 * qJD(2);
	t69 = (qJD(2) * t130 + t109 * t95) * t141;
	t66 = -t69 * t143 - t88 * t95 + (t89 * t131 + t69 * t137) * t114;
	t146 = t66 * t74 * t75;
	t139 = t112 * t84;
	t86 = t104 * t112 - t115 * t133;
	t82 = t86 ^ 2;
	t81 = t82 * t84 + 0.1e1;
	t85 = t83 * t84;
	t124 = -t135 * t116 - t117 * t129;
	t96 = t124 * qJD(2);
	t145 = (-t115 * t82 * t85 + t86 * t139) * t96 / t81 ^ 2;
	t125 = t102 * t109 + t130;
	t70 = t125 * t141;
	t144 = t69 * t70;
	t142 = t124 * t75;
	t140 = t112 * t83;
	t138 = t115 * t86;
	t111 = t109 * t110;
	t99 = t124 ^ 2;
	t97 = t104 * qJD(2);
	t94 = qJD(2) * t100;
	t79 = 0.1e1 / t81;
	t73 = t99 * t75 + 0.1e1;
	t67 = -t70 * t143 - t88 * t102 + (t116 * t89 + t70 * t137) * t114;
	t65 = (-0.2e1 * t125 / t93 ^ 2 * (t100 * t110 * t95 + t111 * t98 * t131) * t108 + (t95 * t134 - t109 * t94 + (t102 * t134 + (0.2e1 * t111 * t116 ^ 2 + t109) * t100) * qJD(2)) * t91) * t107;
	t1 = [0, t65, 0, 0, 0, 0; 0, 0.2e1 * (-t104 * t74 - t67 * t142) / t73 ^ 2 * (-t97 * t142 - t99 * t146) + (-0.2e1 * t67 * t124 * t146 + t96 * t74 + (-t104 * t66 - t67 * t97 - (-(t100 * t144 + t94) * t88 - (-t100 * t65 - t102 * t69 - t70 * t95) * t89) * t124) * t75 + ((-qJD(2) * t70 - t69) * t88 * t116 + (t65 * t88 + (qJD(2) + t144) * t89) * t117) * t114 * t142) / t73, 0, 0, 0, 0; 0, (t84 * t138 - t140) * t97 * t79 - 0.2e1 * (t140 * t145 + (-t84 * t86 * t145 + (-t85 * t138 + t139) * t96 * t79) * t115) * t124, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:32
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (941->56), mult. (2271->131), div. (423->14), fcn. (2956->11), ass. (0->66)
	t146 = sin(qJ(2));
	t147 = cos(qJ(2));
	t144 = sin(pkin(10));
	t173 = cos(pkin(6));
	t161 = t144 * t173;
	t172 = cos(pkin(10));
	t132 = -t146 * t161 + t172 * t147;
	t140 = pkin(11) + qJ(4);
	t136 = sin(t140);
	t137 = cos(t140);
	t145 = sin(pkin(6));
	t166 = t144 * t145;
	t155 = -t132 * t136 + t137 * t166;
	t177 = t155 * qJD(4);
	t158 = t173 * t172;
	t128 = t144 * t146 - t147 * t158;
	t165 = t145 * t147;
	t118 = atan2(-t128, -t165);
	t116 = sin(t118);
	t117 = cos(t118);
	t103 = -t116 * t128 - t117 * t165;
	t100 = 0.1e1 / t103;
	t115 = t132 * t137 + t136 * t166;
	t111 = 0.1e1 / t115;
	t141 = 0.1e1 / t147;
	t101 = 0.1e1 / t103 ^ 2;
	t112 = 0.1e1 / t115 ^ 2;
	t142 = 0.1e1 / t147 ^ 2;
	t130 = t144 * t147 + t146 * t158;
	t123 = t130 * qJD(2);
	t164 = qJD(2) * t146;
	t167 = t142 * t146;
	t162 = t128 * t167;
	t126 = t128 ^ 2;
	t139 = 0.1e1 / t145 ^ 2;
	t121 = t126 * t139 * t142 + 0.1e1;
	t119 = 0.1e1 / t121;
	t138 = 0.1e1 / t145;
	t168 = t119 * t138;
	t95 = (qJD(2) * t162 + t123 * t141) * t168;
	t92 = (-t128 * t95 + t145 * t164) * t117 + (t95 * t165 - t123) * t116;
	t176 = t100 * t101 * t92;
	t110 = t155 ^ 2;
	t106 = t110 * t112 + 0.1e1;
	t154 = -t172 * t146 - t147 * t161;
	t124 = t154 * qJD(2);
	t108 = t115 * qJD(4) + t124 * t136;
	t169 = t112 * t155;
	t109 = t124 * t137 + t177;
	t170 = t109 * t111 * t112;
	t175 = 0.1e1 / t106 ^ 2 * (-t108 * t169 - t110 * t170);
	t156 = t130 * t141 + t162;
	t96 = t156 * t168;
	t174 = t128 * t96;
	t171 = t101 * t154;
	t163 = -0.2e1 * t175;
	t157 = -t111 * t136 - t137 * t169;
	t143 = t141 * t142;
	t127 = t154 ^ 2;
	t125 = t132 * qJD(2);
	t122 = qJD(2) * t128;
	t104 = 0.1e1 / t106;
	t99 = t127 * t101 + 0.1e1;
	t93 = (t145 * t146 - t174) * t117 + (t96 * t165 - t130) * t116;
	t91 = (-0.2e1 * t156 / t121 ^ 2 * (t123 * t128 * t142 + t126 * t143 * t164) * t139 + (t123 * t167 - t122 * t141 + (t130 * t167 + (0.2e1 * t143 * t146 ^ 2 + t141) * t128) * qJD(2)) * t119) * t138;
	t1 = [0, t91, 0, 0, 0, 0; 0, 0.2e1 * (-t100 * t132 - t93 * t171) / t99 ^ 2 * (-t125 * t171 - t127 * t176) + (t124 * t100 + (-t93 * t125 - t132 * t92) * t101 - (0.2e1 * t93 * t176 + (-(-t123 * t96 - t128 * t91 - t130 * t95 + (t95 * t96 + qJD(2)) * t165) * t117 - (t95 * t174 + t122 + (t147 * t91 + (-qJD(2) * t96 - t95) * t146) * t145) * t116) * t101) * t154) / t99, 0, 0, 0, 0; 0, -t157 * t154 * t163 + (t157 * t125 - ((-qJD(4) * t111 + 0.2e1 * t155 * t170) * t137 + (t108 * t137 + (t109 + t177) * t136) * t112) * t154) * t104, 0, t163 - 0.2e1 * (t104 * t108 * t112 - (-t104 * t170 - t112 * t175) * t155) * t155, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:33
	% DurationCPUTime: 1.34s
	% Computational Cost: add. (5351->102), mult. (8196->213), div. (535->12), fcn. (10572->13), ass. (0->102)
	t199 = sin(pkin(10));
	t202 = cos(pkin(10));
	t204 = sin(qJ(2));
	t203 = cos(pkin(6));
	t205 = cos(qJ(2));
	t228 = t203 * t205;
	t187 = -t199 * t204 + t202 * t228;
	t183 = t187 * qJD(2);
	t229 = t203 * t204;
	t188 = t199 * t205 + t202 * t229;
	t197 = pkin(11) + qJ(4);
	t195 = sin(t197);
	t200 = sin(pkin(6));
	t232 = t200 * t202;
	t219 = t195 * t232;
	t196 = cos(t197);
	t225 = qJD(4) * t196;
	t151 = -qJD(4) * t219 + t183 * t195 + t188 * t225;
	t173 = t188 * t195 + t196 * t232;
	t171 = t173 ^ 2;
	t231 = t200 * t204;
	t181 = t195 * t231 - t203 * t196;
	t179 = 0.1e1 / t181 ^ 2;
	t165 = t171 * t179 + 0.1e1;
	t163 = 0.1e1 / t165;
	t182 = t203 * t195 + t196 * t231;
	t226 = qJD(2) * t205;
	t218 = t200 * t226;
	t169 = t182 * qJD(4) + t195 * t218;
	t178 = 0.1e1 / t181;
	t237 = t173 * t179;
	t135 = (-t151 * t178 + t169 * t237) * t163;
	t166 = atan2(-t173, t181);
	t161 = sin(t166);
	t162 = cos(t166);
	t216 = -t161 * t181 - t162 * t173;
	t131 = t216 * t135 - t161 * t151 + t162 * t169;
	t145 = -t161 * t173 + t162 * t181;
	t142 = 0.1e1 / t145;
	t143 = 0.1e1 / t145 ^ 2;
	t251 = t131 * t142 * t143;
	t220 = t199 * t229;
	t190 = t202 * t205 - t220;
	t233 = t199 * t200;
	t177 = t190 * t196 + t195 * t233;
	t198 = sin(pkin(12));
	t189 = t199 * t228 + t202 * t204;
	t201 = cos(pkin(12));
	t234 = t189 * t201;
	t159 = t177 * t198 - t234;
	t185 = t189 * qJD(2);
	t214 = -t190 * t195 + t196 * t233;
	t154 = t214 * qJD(4) - t185 * t196;
	t186 = -qJD(2) * t220 + t202 * t226;
	t150 = t154 * t201 + t186 * t198;
	t235 = t189 * t198;
	t160 = t177 * t201 + t235;
	t156 = 0.1e1 / t160;
	t157 = 0.1e1 / t160 ^ 2;
	t245 = t150 * t156 * t157;
	t250 = 0.2e1 * t159 * t245;
	t249 = -0.2e1 * t214 * t251;
	t230 = t200 * t205;
	t212 = -t178 * t187 + t230 * t237;
	t248 = t195 * t212;
	t238 = t169 * t178 * t179;
	t247 = -0.2e1 * (t151 * t237 - t171 * t238) / t165 ^ 2;
	t246 = t143 * t214;
	t153 = t177 * qJD(4) - t185 * t195;
	t244 = t153 * t143;
	t243 = t156 * t198;
	t242 = t157 * t159;
	t241 = t159 * t201;
	t240 = t161 * t214;
	t239 = t162 * t214;
	t236 = t189 * t195;
	t227 = qJD(2) * t204;
	t172 = t214 ^ 2;
	t141 = t172 * t143 + 0.1e1;
	t224 = 0.2e1 * (-t172 * t251 - t214 * t244) / t141 ^ 2;
	t155 = t159 ^ 2;
	t148 = t155 * t157 + 0.1e1;
	t149 = t154 * t198 - t186 * t201;
	t223 = 0.2e1 * (t149 * t242 - t155 * t245) / t148 ^ 2;
	t217 = -0.2e1 * t173 * t238;
	t175 = t188 * t196 - t219;
	t215 = -t175 * t178 + t182 * t237;
	t213 = qJD(4) * t236 - t186 * t196;
	t184 = t188 * qJD(2);
	t170 = -t181 * qJD(4) + t196 * t218;
	t168 = t190 * t198 - t196 * t234;
	t167 = -t190 * t201 - t196 * t235;
	t152 = -t173 * qJD(4) + t183 * t196;
	t146 = 0.1e1 / t148;
	t138 = 0.1e1 / t141;
	t137 = t163 * t248;
	t136 = t215 * t163;
	t133 = (-t161 * t187 + t162 * t230) * t195 + t216 * t137;
	t132 = t216 * t136 - t161 * t175 + t162 * t182;
	t130 = t215 * t247 + (t182 * t217 - t152 * t178 + (t151 * t182 + t169 * t175 + t170 * t173) * t179) * t163;
	t128 = t247 * t248 + (t212 * t225 + (t217 * t230 + t178 * t184 + (t169 * t187 + (t151 * t205 - t173 * t227) * t200) * t179) * t195) * t163;
	t1 = [0, t128, 0, t130, 0, 0; 0, (-t133 * t246 + t142 * t236) * t224 + ((-t186 * t195 - t189 * t225) * t142 + (-t244 + t249) * t133 + (t236 * t131 + (-t128 * t173 - t137 * t151 + (-t195 * t227 + t205 * t225) * t200 + (-t137 * t181 - t187 * t195) * t135) * t239 + (-t187 * t225 - t128 * t181 - t137 * t169 + t184 * t195 + (t137 * t173 - t195 * t230) * t135) * t240) * t143) * t138, 0, (-t132 * t246 - t142 * t177) * t224 + (t132 * t249 + t154 * t142 + (-t177 * t131 - t132 * t153 + (-t130 * t173 - t136 * t151 + t170 + (-t136 * t181 - t175) * t135) * t239 + (-t130 * t181 - t136 * t169 - t152 + (t136 * t173 - t182) * t135) * t240) * t143) * t138, 0, 0; 0, (-t156 * t167 + t168 * t242) * t223 + ((t185 * t201 + t198 * t213) * t156 + t168 * t250 + (-t167 * t150 - (-t185 * t198 + t201 * t213) * t159 - t168 * t149) * t157) * t146, 0, -(-t157 * t241 + t243) * t214 * t223 + (t214 * t201 * t250 - t153 * t243 + (t153 * t241 - (t149 * t201 + t150 * t198) * t214) * t157) * t146, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:32
	% EndTime: 2019-10-09 21:35:33
	% DurationCPUTime: 1.45s
	% Computational Cost: add. (6115->111), mult. (9085->227), div. (559->12), fcn. (11668->13), ass. (0->106)
	t230 = sin(pkin(10));
	t232 = cos(pkin(10));
	t234 = sin(qJ(2));
	t233 = cos(pkin(6));
	t235 = cos(qJ(2));
	t261 = t233 * t235;
	t216 = -t230 * t234 + t232 * t261;
	t212 = t216 * qJD(2);
	t262 = t233 * t234;
	t217 = t230 * t235 + t232 * t262;
	t229 = pkin(11) + qJ(4);
	t225 = sin(t229);
	t231 = sin(pkin(6));
	t265 = t231 * t232;
	t251 = t225 * t265;
	t227 = cos(t229);
	t258 = qJD(4) * t227;
	t185 = -qJD(4) * t251 + t212 * t225 + t217 * t258;
	t201 = t217 * t225 + t227 * t265;
	t199 = t201 ^ 2;
	t264 = t231 * t234;
	t210 = t225 * t264 - t233 * t227;
	t208 = 0.1e1 / t210 ^ 2;
	t193 = t199 * t208 + 0.1e1;
	t191 = 0.1e1 / t193;
	t211 = t233 * t225 + t227 * t264;
	t259 = qJD(2) * t235;
	t250 = t231 * t259;
	t197 = t211 * qJD(4) + t225 * t250;
	t207 = 0.1e1 / t210;
	t269 = t201 * t208;
	t163 = (-t185 * t207 + t197 * t269) * t191;
	t194 = atan2(-t201, t210);
	t189 = sin(t194);
	t190 = cos(t194);
	t247 = -t189 * t210 - t190 * t201;
	t159 = t247 * t163 - t189 * t185 + t190 * t197;
	t173 = -t189 * t201 + t190 * t210;
	t170 = 0.1e1 / t173;
	t171 = 0.1e1 / t173 ^ 2;
	t283 = t159 * t170 * t171;
	t252 = t230 * t262;
	t219 = t232 * t235 - t252;
	t266 = t230 * t231;
	t244 = -t219 * t225 + t227 * t266;
	t282 = -0.2e1 * t244 * t283;
	t263 = t231 * t235;
	t243 = -t207 * t216 + t263 * t269;
	t281 = t225 * t243;
	t270 = t197 * t207 * t208;
	t280 = -0.2e1 * (t185 * t269 - t199 * t270) / t193 ^ 2;
	t205 = t219 * t227 + t225 * t266;
	t218 = t230 * t261 + t232 * t234;
	t228 = pkin(12) + qJ(6);
	t224 = sin(t228);
	t226 = cos(t228);
	t184 = t205 * t226 + t218 * t224;
	t180 = 0.1e1 / t184;
	t181 = 0.1e1 / t184 ^ 2;
	t214 = t218 * qJD(2);
	t188 = t244 * qJD(4) - t214 * t227;
	t215 = -qJD(2) * t252 + t232 * t259;
	t174 = t184 * qJD(6) + t188 * t224 - t215 * t226;
	t183 = t205 * t224 - t218 * t226;
	t179 = t183 ^ 2;
	t178 = t179 * t181 + 0.1e1;
	t275 = t181 * t183;
	t257 = qJD(6) * t183;
	t175 = t188 * t226 + t215 * t224 - t257;
	t277 = t175 * t180 * t181;
	t279 = (t174 * t275 - t179 * t277) / t178 ^ 2;
	t278 = t171 * t244;
	t276 = t180 * t224;
	t274 = t183 * t226;
	t187 = t205 * qJD(4) - t214 * t225;
	t273 = t187 * t171;
	t272 = t189 * t244;
	t271 = t190 * t244;
	t268 = t218 * t225;
	t267 = t218 * t227;
	t260 = qJD(2) * t234;
	t200 = t244 ^ 2;
	t169 = t200 * t171 + 0.1e1;
	t256 = 0.2e1 * (-t200 * t283 - t244 * t273) / t169 ^ 2;
	t255 = -0.2e1 * t279;
	t253 = t183 * t277;
	t249 = -0.2e1 * t201 * t270;
	t248 = qJD(6) * t267 - t214;
	t246 = t181 * t274 - t276;
	t203 = t217 * t227 - t251;
	t245 = -t203 * t207 + t211 * t269;
	t242 = qJD(4) * t268 + qJD(6) * t219 - t215 * t227;
	t213 = t217 * qJD(2);
	t198 = -t210 * qJD(4) + t227 * t250;
	t196 = t219 * t224 - t226 * t267;
	t195 = -t219 * t226 - t224 * t267;
	t186 = -t201 * qJD(4) + t212 * t227;
	t176 = 0.1e1 / t178;
	t166 = 0.1e1 / t169;
	t165 = t191 * t281;
	t164 = t245 * t191;
	t161 = (-t189 * t216 + t190 * t263) * t225 + t247 * t165;
	t160 = t247 * t164 - t189 * t203 + t190 * t211;
	t158 = t245 * t280 + (t211 * t249 - t186 * t207 + (t185 * t211 + t197 * t203 + t198 * t201) * t208) * t191;
	t156 = t280 * t281 + (t243 * t258 + (t249 * t263 + t207 * t213 + (t197 * t216 + (t185 * t235 - t201 * t260) * t231) * t208) * t225) * t191;
	t1 = [0, t156, 0, t158, 0, 0; 0, (-t161 * t278 + t170 * t268) * t256 + ((-t215 * t225 - t218 * t258) * t170 + (-t273 + t282) * t161 + (t268 * t159 + (-t156 * t201 - t165 * t185 + (-t225 * t260 + t235 * t258) * t231 + (-t165 * t210 - t216 * t225) * t163) * t271 + (-t216 * t258 - t156 * t210 - t165 * t197 + t213 * t225 + (t165 * t201 - t225 * t263) * t163) * t272) * t171) * t166, 0, (-t160 * t278 - t170 * t205) * t256 + (t160 * t282 + t188 * t170 + (-t205 * t159 - t160 * t187 + (-t158 * t201 - t164 * t185 + t198 + (-t164 * t210 - t203) * t163) * t271 + (-t158 * t210 - t164 * t197 - t186 + (t164 * t201 - t211) * t163) * t272) * t171) * t166, 0, 0; 0, 0.2e1 * (-t180 * t195 + t196 * t275) * t279 + (0.2e1 * t196 * t253 - t248 * t180 * t226 + t242 * t276 + (-t248 * t183 * t224 - t196 * t174 - t195 * t175 - t242 * t274) * t181) * t176, 0, -t246 * t244 * t255 + (t246 * t187 - ((-qJD(6) * t180 - 0.2e1 * t253) * t226 + (t174 * t226 + (t175 - t257) * t224) * t181) * t244) * t176, 0, t255 + 0.2e1 * (t174 * t181 * t176 + (-t176 * t277 - t181 * t279) * t183) * t183;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end