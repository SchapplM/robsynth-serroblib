% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRPP3
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
%   Wie in S5PRRPP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRRPP3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPP3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:14:24
	% EndTime: 2019-12-05 16:14:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:14:24
	% EndTime: 2019-12-05 16:14:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:14:24
	% EndTime: 2019-12-05 16:14:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:14:24
	% EndTime: 2019-12-05 16:14:25
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (429->41), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->59)
	t97 = sin(qJ(2));
	t90 = t97 ^ 2;
	t99 = cos(qJ(2));
	t123 = t90 / t99 ^ 2;
	t111 = 0.1e1 + t123;
	t94 = sin(pkin(7));
	t135 = t111 * t94;
	t134 = qJD(2) * t97;
	t95 = cos(pkin(7));
	t120 = t95 * t99;
	t96 = sin(qJ(3));
	t98 = cos(qJ(3));
	t77 = t98 * t120 + t94 * t96;
	t121 = t94 * t97;
	t80 = atan2(-t121, -t99);
	t78 = sin(t80);
	t79 = cos(t80);
	t63 = -t121 * t78 - t79 * t99;
	t60 = 0.1e1 / t63;
	t73 = 0.1e1 / t77;
	t124 = t78 * t99;
	t108 = -t124 * t94 + t79 * t97;
	t115 = t79 * t121;
	t109 = -t115 + t124;
	t87 = t94 ^ 2;
	t83 = t123 * t87 + 0.1e1;
	t81 = 0.1e1 / t83;
	t67 = t81 * t135;
	t54 = t109 * t67 + t108;
	t132 = 0.2e1 * t54;
	t61 = 0.1e1 / t63 ^ 2;
	t74 = 0.1e1 / t77 ^ 2;
	t118 = qJD(2) * t99;
	t127 = t61 * t97;
	t64 = qJD(2) * t67;
	t53 = qJD(2) * t108 + t109 * t64;
	t130 = t53 * t60 * t61;
	t88 = t95 ^ 2;
	t59 = t61 * t88 * t90 + 0.1e1;
	t131 = (t118 * t127 - t130 * t90) * t88 / t59 ^ 2;
	t76 = t120 * t96 - t94 * t98;
	t125 = t74 * t76;
	t113 = t95 * t134;
	t117 = qJD(3) * t76;
	t70 = -t113 * t98 - t117;
	t126 = t70 * t73 * t74;
	t72 = t76 ^ 2;
	t68 = t72 * t74 + 0.1e1;
	t69 = -qJD(3) * t77 + t96 * t113;
	t129 = (-t125 * t69 - t126 * t72) / t68 ^ 2;
	t57 = 0.1e1 / t59;
	t128 = t57 * t61;
	t119 = t67 - t94;
	t116 = -0.2e1 * t129;
	t112 = t67 * t94 - 0.1e1;
	t110 = t125 * t98 - t73 * t96;
	t65 = 0.1e1 / t68;
	t56 = 0.2e1 * (t81 - t111 / t83 ^ 2 * t87) / t99 * t134 * t135;
	t1 = [0, t56, 0, 0, 0; 0, ((-0.2e1 * t60 * t131 + (-qJD(2) * t54 - t53) * t128) * t99 + (t61 * t131 * t132 + (t56 * t61 * t115 + t130 * t132 - qJD(2) * t60 - (-qJD(2) * t119 + t112 * t64) * t78 * t127) * t57 - (t56 * t78 + (-qJD(2) * t112 + t119 * t64) * t79) * t99 * t128) * t97) * t95, 0, 0, 0; 0, (t110 * t97 * t116 + (t110 * t118 + ((-qJD(3) * t73 - 0.2e1 * t126 * t76) * t98 + (-t69 * t98 + (t70 - t117) * t96) * t74) * t97) * t65) * t95, t116 + 0.2e1 * (-t65 * t69 * t74 + (-t126 * t65 - t129 * t74) * t76) * t76, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:14:24
	% EndTime: 2019-12-05 16:14:25
	% DurationCPUTime: 0.80s
	% Computational Cost: add. (1276->95), mult. (4381->222), div. (832->14), fcn. (5598->11), ass. (0->96)
	t146 = sin(pkin(7));
	t148 = cos(pkin(7));
	t151 = cos(qJ(3));
	t149 = sin(qJ(3));
	t152 = cos(qJ(2));
	t182 = t149 * t152;
	t129 = t146 * t182 + t148 * t151;
	t150 = sin(qJ(2));
	t181 = t150 * t149;
	t123 = atan2(-t129, t181);
	t119 = sin(t123);
	t120 = cos(t123);
	t127 = t129 ^ 2;
	t140 = 0.1e1 / t149 ^ 2;
	t143 = 0.1e1 / t150 ^ 2;
	t186 = t140 * t143;
	t124 = t127 * t186 + 0.1e1;
	t121 = 0.1e1 / t124;
	t139 = 0.1e1 / t149;
	t171 = t139 * t143 * t152;
	t163 = t129 * t171 + t146;
	t99 = t163 * t121;
	t197 = t146 - t99;
	t204 = t119 * t150 * t197 + t120 * t152;
	t190 = t119 * t129;
	t103 = t120 * t181 - t190;
	t101 = 0.1e1 / t103 ^ 2;
	t175 = qJD(3) * t151;
	t165 = t152 * t175;
	t178 = qJD(2) * t150;
	t115 = -t146 * t165 + (qJD(3) * t148 + t146 * t178) * t149;
	t177 = qJD(2) * t152;
	t167 = t149 * t177;
	t161 = t150 * t175 + t167;
	t142 = 0.1e1 / t150;
	t172 = t129 * t186;
	t93 = (t115 * t142 * t139 + t161 * t172) * t121;
	t198 = t129 * t93;
	t89 = (-t93 * t181 + t115) * t119 + (t161 - t198) * t120;
	t203 = t89 * t101;
	t179 = t151 * t152;
	t184 = t148 * t149;
	t131 = t146 * t179 - t184;
	t185 = t140 * t151;
	t162 = t129 * t185 - t131 * t139;
	t202 = t142 * t162;
	t168 = t151 * t177;
	t176 = qJD(3) * t149;
	t201 = t150 * t176 - t168;
	t141 = t139 * t140;
	t144 = t142 * t143;
	t200 = -0.2e1 / t124 ^ 2 * (-t115 * t172 + (-t140 * t144 * t177 - t141 * t143 * t175) * t127);
	t100 = 0.1e1 / t103;
	t133 = t146 * t149 + t148 * t179;
	t145 = sin(pkin(8));
	t147 = cos(pkin(8));
	t183 = t148 * t150;
	t114 = t133 * t147 + t145 * t183;
	t110 = 0.1e1 / t114;
	t111 = 0.1e1 / t114 ^ 2;
	t199 = t100 * t203;
	t132 = -t146 * t151 + t148 * t182;
	t196 = t101 * t132;
	t169 = t151 * t178;
	t118 = -t132 * qJD(3) - t148 * t169;
	t170 = t148 * t177;
	t108 = t118 * t147 + t145 * t170;
	t112 = t110 * t111;
	t195 = t108 * t112;
	t194 = t110 * t145;
	t113 = t133 * t145 - t147 * t183;
	t193 = t111 * t113;
	t180 = t150 * t151;
	t126 = (t145 * t152 - t147 * t180) * t148;
	t192 = t113 * t126;
	t191 = t113 * t147;
	t188 = t120 * t129;
	t117 = -t146 * t176 - t148 * t165 + t178 * t184;
	t128 = t132 ^ 2;
	t98 = t128 * t101 + 0.1e1;
	t174 = 0.2e1 * (-t117 * t196 - t128 * t199) / t98 ^ 2;
	t109 = t113 ^ 2;
	t106 = t109 * t111 + 0.1e1;
	t107 = t118 * t145 - t147 * t170;
	t173 = 0.2e1 / t106 ^ 2 * (t107 * t193 - t109 * t195);
	t164 = 0.2e1 * t132 * t199;
	t125 = (-t145 * t180 - t147 * t152) * t148;
	t116 = qJD(3) * t129 + t146 * t169;
	t104 = 0.1e1 / t106;
	t96 = 0.1e1 / t98;
	t95 = t121 * t202;
	t91 = t204 * t149 - t99 * t188;
	t90 = (-t129 * t95 + t180) * t120 + (-t95 * t181 - t131) * t119;
	t88 = t163 * t200 + (-t115 * t171 + (-t165 * t186 + (-0.2e1 * t144 * t152 ^ 2 - t142) * t139 * qJD(2)) * t129) * t121;
	t86 = t200 * t202 + (-t162 * t143 * t177 + (-t115 * t185 + t116 * t139 + (t131 * t185 + (-0.2e1 * t141 * t151 ^ 2 - t139) * t129) * qJD(3)) * t142) * t121;
	t1 = [0, t88, t86, 0, 0; 0, t91 * t174 * t196 + (t91 * t164 + (-(-t88 * t188 + (t115 * t120 + t190 * t93) * t99) * t132 + t91 * t117) * t101 + (-t100 * t183 - t204 * t196) * t175) * t96 + (-((t197 * t93 - qJD(2)) * t150 * t120 + (-t150 * t88 + (t197 * qJD(2) - t93) * t152) * t119) * t96 * t196 + (t150 * t96 * t203 + (t150 * t174 - t96 * t177) * t100) * t148) * t149, (-t100 * t133 + t90 * t196) * t174 + (t90 * t164 + t118 * t100 + (t90 * t117 - t133 * t89) * t101 + (-(t168 + t115 * t95 - t129 * t86 - t131 * t93 + (-t93 * t95 - qJD(3)) * t181) * t120 - (t116 + (-t167 + t198) * t95 + (-t149 * t86 + (-qJD(3) * t95 - t93) * t151) * t150) * t119) * t196) * t96, 0, 0; 0, (-t110 * t125 + t111 * t192) * t173 + (-t126 * t107 * t111 + (-t111 * t125 + 0.2e1 * t112 * t192) * t108 + ((t201 * t145 + t147 * t178) * t110 - (-t145 * t178 + t201 * t147) * t193) * t148) * t104, (-t111 * t191 + t194) * t132 * t173 + (-0.2e1 * t132 * t191 * t195 + t117 * t194 + (-t117 * t191 + (t107 * t147 + t108 * t145) * t132) * t111) * t104, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:14:24
	% EndTime: 2019-12-05 16:14:25
	% DurationCPUTime: 0.80s
	% Computational Cost: add. (2314->93), mult. (7265->210), div. (523->12), fcn. (9299->11), ass. (0->100)
	t162 = sin(pkin(8));
	t168 = cos(qJ(2));
	t167 = cos(qJ(3));
	t220 = sin(pkin(7));
	t188 = t220 * t167;
	t182 = t168 * t188;
	t164 = cos(pkin(7));
	t165 = sin(qJ(3));
	t206 = t164 * t165;
	t175 = t182 - t206;
	t163 = cos(pkin(8));
	t166 = sin(qJ(2));
	t190 = t166 * t220;
	t184 = t163 * t190;
	t138 = t175 * t162 - t184;
	t133 = t138 ^ 2;
	t201 = t168 * t163;
	t203 = t166 * t167;
	t155 = t162 * t203 + t201;
	t151 = 0.1e1 / t155 ^ 2;
	t127 = t133 * t151 + 0.1e1;
	t189 = t220 * t165;
	t154 = -t164 * t167 - t168 * t189;
	t181 = qJD(2) * t190;
	t183 = t220 * t201;
	t129 = (t154 * qJD(3) - t167 * t181) * t162 - qJD(2) * t183;
	t202 = t167 * t168;
	t204 = t166 * t163;
	t158 = t162 * t202 - t204;
	t197 = qJD(3) * t166;
	t146 = -t165 * t162 * t197 + t158 * qJD(2);
	t150 = 0.1e1 / t155;
	t209 = t146 * t150 * t151;
	t211 = t138 * t151;
	t221 = (t129 * t211 - t133 * t209) / t127 ^ 2;
	t128 = atan2(-t138, t155);
	t120 = sin(t128);
	t121 = cos(t128);
	t179 = t120 * t138 - t121 * t155;
	t116 = 0.1e1 / t179;
	t157 = t164 * t202 + t189;
	t207 = t162 * t164;
	t141 = t157 * t163 + t166 * t207;
	t135 = 0.1e1 / t141;
	t117 = 0.1e1 / t179 ^ 2;
	t136 = 0.1e1 / t141 ^ 2;
	t124 = 0.1e1 / t127;
	t108 = (-t129 * t150 + t146 * t211) * t124;
	t180 = -t120 * t155 - t121 * t138;
	t105 = t180 * t108 - t120 * t129 + t121 * t146;
	t219 = t105 * t116 * t117;
	t156 = -t168 * t206 + t188;
	t153 = t156 ^ 2;
	t126 = t153 * t136 + 0.1e1;
	t200 = qJD(2) * t166;
	t143 = -t157 * qJD(3) + t200 * t206;
	t212 = t136 * t156;
	t193 = t143 * t212;
	t198 = qJD(3) * t165;
	t191 = t168 * t198;
	t144 = -qJD(3) * t188 + (t167 * t200 + t191) * t164;
	t199 = qJD(2) * t168;
	t131 = -t144 * t163 + t199 * t207;
	t137 = t135 * t136;
	t213 = t131 * t137;
	t218 = (-t153 * t213 + t193) / t126 ^ 2;
	t140 = t157 * t162 - t164 * t204;
	t217 = t117 * t140;
	t216 = t120 * t140;
	t215 = t121 * t140;
	t192 = t163 * t199;
	t130 = -t144 * t162 - t164 * t192;
	t214 = t130 * t117;
	t210 = t138 * t158;
	t208 = t153 * t163;
	t205 = t165 * t166;
	t134 = t140 ^ 2;
	t113 = t134 * t117 + 0.1e1;
	t196 = 0.2e1 * (t134 * t219 + t140 * t214) / t113 ^ 2;
	t195 = 0.2e1 * t218;
	t194 = t138 * t205;
	t187 = -0.2e1 * t140 * t219;
	t185 = t162 * t190;
	t147 = t167 * t185 + t183;
	t178 = t147 * t150 + t151 * t210;
	t177 = t150 * t154 + t151 * t194;
	t176 = t165 * t199 + t167 * t197;
	t148 = t155 * t164;
	t145 = -t155 * qJD(2) - t162 * t191;
	t142 = -t175 * qJD(3) + t165 * t181;
	t132 = -t185 * t198 + (t162 * t182 - t184) * qJD(2);
	t122 = 0.1e1 / t126;
	t111 = 0.1e1 / t113;
	t110 = t177 * t162 * t124;
	t109 = t178 * t124;
	t107 = (-t120 * t154 - t121 * t205) * t162 - t180 * t110;
	t106 = t180 * t109 + t120 * t147 + t121 * t158;
	t104 = -0.2e1 * t178 * t221 + (-0.2e1 * t209 * t210 + t132 * t150 + (t129 * t158 + t138 * t145 - t146 * t147) * t151) * t124;
	t102 = (0.2e1 * t177 * t221 + (0.2e1 * t194 * t209 - t142 * t150 + (-t129 * t205 - t176 * t138 + t146 * t154) * t151) * t124) * t162;
	t1 = [0, t104, t102, 0, 0; 0, (t106 * t217 - t116 * t148) * t196 + (t106 * t187 + t146 * t116 * t164 + (t148 * t105 - t106 * t130 - (-t104 * t138 - t109 * t129 + t145 + (-t109 * t155 + t147) * t108) * t215 - (-t104 * t155 - t109 * t146 + t132 + (t109 * t138 - t158) * t108) * t216) * t117) * t111, (t116 * t156 * t162 + t107 * t217) * t196 + (-(t180 * t102 - (t179 * t108 - t120 * t146 - t121 * t129) * t110) * t217 + (t187 - t214) * t107 + (-t143 * t116 + (-t156 * t105 - (t108 * t205 - t142) * t216 - (-t108 * t154 - t176) * t215) * t117) * t162) * t111, 0, 0; 0, ((t195 * t212 + (-t136 * t143 + 0.2e1 * t156 * t213) * t122) * (t162 * t168 - t163 * t203) - 0.2e1 * t135 * t205 * t218 + ((t167 * t156 * t192 + (-t165 * t131 - (-qJD(2) * t162 + t163 * t198) * t156) * t166) * t136 + t176 * t135) * t122) * t164, (t135 * t157 + t136 * t208) * t195 + (-0.2e1 * t163 * t193 + t135 * t144 + (t136 * t157 + 0.2e1 * t137 * t208) * t131) * t122, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end