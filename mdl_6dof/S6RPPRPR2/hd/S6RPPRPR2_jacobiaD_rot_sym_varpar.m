% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR2
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
%   Wie in S6RPPRPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:36
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:01
	% DurationCPUTime: 0.79s
	% Computational Cost: add. (2561->72), mult. (1839->155), div. (436->14), fcn. (2165->7), ass. (0->77)
	t99 = qJ(1) + pkin(9);
	t97 = cos(t99);
	t128 = qJD(1) * t97;
	t95 = sin(t99);
	t150 = 0.2e1 * t95;
	t84 = t95 ^ 2;
	t85 = 0.1e1 / t95;
	t93 = t97 ^ 2;
	t148 = (0.1e1 + 0.1e1 / t84 * t93) * t85 * t128;
	t98 = pkin(10) + qJ(4);
	t94 = sin(t98);
	t130 = t95 * t94;
	t96 = cos(t98);
	t75 = atan2(-t130, -t96);
	t73 = sin(t75);
	t120 = t73 * t130;
	t74 = cos(t75);
	t69 = -t74 * t96 - t120;
	t66 = 0.1e1 / t69;
	t89 = 0.1e1 / t96;
	t147 = -0.2e1 * t94;
	t67 = 0.1e1 / t69 ^ 2;
	t83 = t94 ^ 2;
	t90 = 0.1e1 / t96 ^ 2;
	t135 = t83 * t90;
	t80 = t84 * t135 + 0.1e1;
	t76 = 0.1e1 / t80;
	t146 = t76 - 0.1e1;
	t118 = t94 * t128;
	t127 = qJD(4) * t95;
	t137 = t74 * t94;
	t126 = qJD(4) * t96;
	t62 = (-(-t95 * t126 - t118) * t89 + t127 * t135) * t76;
	t58 = (-t62 * t95 + qJD(4)) * t137 + (-t118 + (t62 - t127) * t96) * t73;
	t145 = t58 * t66 * t67;
	t144 = t62 * t73;
	t143 = t62 * t94;
	t142 = t67 * t94;
	t141 = t67 * t97;
	t132 = t89 * t94;
	t82 = t94 * t83;
	t88 = t96 ^ 2;
	t106 = qJD(4) * (t82 * t89 / t88 + t132);
	t111 = t83 * t95 * t128;
	t140 = (t84 * t106 + t90 * t111) / t80 ^ 2;
	t117 = 0.1e1 + t135;
	t72 = t117 * t95 * t76;
	t139 = t72 * t95;
	t138 = t73 * t96;
	t136 = t83 * t89;
	t134 = t83 * t93;
	t86 = 0.1e1 / t95 ^ 2;
	t133 = t86 * t93;
	t131 = t94 * t97;
	t129 = qJD(1) * t95;
	t125 = qJD(4) * t97;
	t110 = t93 * t94 * t126;
	t65 = t67 * t134 + 0.1e1;
	t124 = 0.2e1 * (-t134 * t145 + (t110 - t111) * t67) / t65 ^ 2;
	t123 = 0.2e1 * t145;
	t81 = t88 * t133 + 0.1e1;
	t122 = 0.2e1 * (-t86 * t110 - t88 * t148) / t81 ^ 2;
	t121 = t67 * t131;
	t119 = t76 * t136;
	t116 = 0.1e1 + t133;
	t115 = t94 * t124;
	t114 = t140 * t147;
	t113 = t140 * t150;
	t112 = t95 * t119;
	t109 = t117 * t97;
	t108 = t116 * t94;
	t78 = 0.1e1 / t81;
	t63 = 0.1e1 / t65;
	t61 = (t146 * t94 * t73 - t74 * t112) * t97;
	t60 = -t95 * t138 + t137 + (-t74 * t130 + t138) * t72;
	t59 = -t117 * t113 + (qJD(1) * t109 + t106 * t150) * t76;
	t1 = [t97 * t89 * t114 + (qJD(4) * t109 - t129 * t132) * t76, 0, 0, t59, 0, 0; (t66 * t115 + (-t66 * t126 + (qJD(1) * t61 + t58) * t142) * t63) * t95 + (t67 * t115 * t61 + (-((t62 * t112 + t146 * t126 + t114) * t73 + (t113 * t136 - t143 + (t143 + (-t82 * t90 + t147) * t127) * t76) * t74) * t121 + (t94 * t123 - t67 * t126) * t61 + (-t66 + ((-t84 + t93) * t74 * t119 + t146 * t120) * t67) * t94 * qJD(1)) * t63) * t97, 0, 0, (t60 * t142 - t66 * t96) * t97 * t124 + ((-t66 * t129 + (-qJD(4) * t60 - t58) * t141) * t96 + (-t66 * t125 - (-t59 * t74 * t95 + t73 * t127 + t139 * t144 - t144 + (-qJD(4) * t73 - t128 * t74) * t72) * t121 + (t97 * t123 + t67 * t129) * t60 - ((t59 - t128) * t73 + ((0.1e1 - t139) * qJD(4) + (t72 - t95) * t62) * t74) * t96 * t141) * t94) * t63, 0, 0; t116 * t96 * t122 + (qJD(4) * t108 + 0.2e1 * t96 * t148) * t78, 0, 0, t85 * t122 * t131 + (-t85 * t96 * t125 + qJD(1) * t108) * t78, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:01
	% EndTime: 2019-10-09 23:36:02
	% DurationCPUTime: 1.06s
	% Computational Cost: add. (3055->92), mult. (2519->204), div. (480->12), fcn. (2968->9), ass. (0->94)
	t142 = pkin(10) + qJ(4);
	t138 = sin(t142);
	t132 = 0.1e1 / t138 ^ 2;
	t140 = cos(t142);
	t136 = t140 ^ 2;
	t192 = t132 * t136;
	t207 = t140 * t192;
	t143 = qJ(1) + pkin(9);
	t139 = sin(t143);
	t165 = 0.1e1 + t192;
	t206 = t139 * t165;
	t144 = sin(qJ(6));
	t145 = cos(qJ(6));
	t162 = qJD(6) * t138 + qJD(1);
	t177 = qJD(4) * t140;
	t205 = t162 * t144 - t145 * t177;
	t204 = t144 * t177 + t162 * t145;
	t189 = t139 * t140;
	t126 = atan2(-t189, t138);
	t125 = cos(t126);
	t124 = sin(t126);
	t171 = t124 * t189;
	t112 = t125 * t138 - t171;
	t108 = 0.1e1 / t112;
	t141 = cos(t143);
	t185 = t141 * t144;
	t187 = t139 * t145;
	t121 = t138 * t185 + t187;
	t117 = 0.1e1 / t121;
	t131 = 0.1e1 / t138;
	t109 = 0.1e1 / t112 ^ 2;
	t118 = 0.1e1 / t121 ^ 2;
	t134 = t139 ^ 2;
	t129 = t134 * t192 + 0.1e1;
	t127 = 0.1e1 / t129;
	t203 = t127 - 0.1e1;
	t180 = qJD(1) * t141;
	t169 = t140 * t180;
	t178 = qJD(4) * t139;
	t101 = ((t138 * t178 - t169) * t131 + t178 * t192) * t127;
	t193 = t125 * t140;
	t96 = (-t101 * t139 + qJD(4)) * t193 + (-t169 + (-t101 + t178) * t138) * t124;
	t202 = t108 * t109 * t96;
	t161 = qJD(1) * t138 + qJD(6);
	t156 = t161 * t145;
	t105 = t139 * t156 + t205 * t141;
	t184 = t141 * t145;
	t188 = t139 * t144;
	t120 = -t138 * t184 + t188;
	t116 = t120 ^ 2;
	t115 = t116 * t118 + 0.1e1;
	t195 = t118 * t120;
	t157 = t161 * t144;
	t106 = -t139 * t157 + t204 * t141;
	t199 = t106 * t117 * t118;
	t201 = 0.1e1 / t115 ^ 2 * (t105 * t195 - t116 * t199);
	t200 = t101 * t140;
	t154 = qJD(4) * (-t140 - t207) * t131;
	t190 = t136 * t139;
	t159 = t180 * t190;
	t198 = (t132 * t159 + t134 * t154) / t129 ^ 2;
	t197 = t109 * t140;
	t196 = t109 * t141;
	t194 = t124 * t138;
	t137 = t141 ^ 2;
	t191 = t136 * t137;
	t186 = t140 * t141;
	t111 = t127 * t206;
	t183 = -t111 + t139;
	t182 = qJD(1) * t139;
	t181 = qJD(1) * t140;
	t179 = qJD(4) * t138;
	t104 = t109 * t191 + 0.1e1;
	t176 = 0.2e1 / t104 ^ 2 * (-t191 * t202 + (-t137 * t138 * t177 - t159) * t109);
	t175 = 0.2e1 * t202;
	t174 = 0.2e1 * t201;
	t173 = -0.2e1 * t198;
	t172 = t140 * t198;
	t170 = t131 * t190;
	t166 = t111 * t139 - 0.1e1;
	t164 = t140 * t176;
	t163 = 0.2e1 * t120 * t199;
	t160 = t127 * t170;
	t158 = t165 * t141;
	t155 = t117 * t145 + t144 * t195;
	t153 = t155 * t141;
	t123 = -t138 * t188 + t184;
	t122 = t138 * t187 + t185;
	t113 = 0.1e1 / t115;
	t102 = 0.1e1 / t104;
	t100 = (t203 * t140 * t124 + t125 * t160) * t141;
	t98 = t139 * t194 + t193 + (-t125 * t189 - t194) * t111;
	t97 = t173 * t206 + (qJD(1) * t158 + 0.2e1 * t139 * t154) * t127;
	t1 = [0.2e1 * t131 * t141 * t172 + (t131 * t139 * t181 + qJD(4) * t158) * t127, 0, 0, t97, 0, 0; (t108 * t164 + (t108 * t179 + (qJD(1) * t100 + t96) * t197) * t102) * t139 + (t109 * t164 * t100 + (-((-t101 * t160 - t203 * t179 - 0.2e1 * t172) * t124 + (t170 * t173 - t200 + (t200 + (-0.2e1 * t140 - t207) * t178) * t127) * t125) * t109 * t186 + (t109 * t179 + t140 * t175) * t100 + (-t108 + ((t134 - t137) * t136 * t131 * t127 * t125 + t203 * t171) * t109) * t181) * t102) * t141, 0, 0, (t108 * t138 + t98 * t197) * t141 * t176 + ((t108 * t182 + (qJD(4) * t98 + t96) * t196) * t138 + ((-qJD(4) * t108 + t98 * t175) * t141 + (t98 * t182 + (-(-t111 * t180 - t139 * t97) * t125 - (t183 * qJD(4) + t166 * t101) * t124) * t186) * t109 - ((-t97 + t180) * t124 + (t166 * qJD(4) + t183 * t101) * t125) * t138 * t196) * t140) * t102, 0, 0; (-t117 * t122 + t123 * t195) * t174 + (t123 * t163 + (-t123 * t105 - t122 * t106 + (t204 * t139 + t141 * t157) * t120) * t118 + (-t205 * t139 + t141 * t156) * t117) * t113, 0, 0, t140 * t153 * t174 + (t153 * t179 + (t155 * t182 + ((qJD(6) * t117 + t163) * t144 + (-t105 * t144 + (-qJD(6) * t120 + t106) * t145) * t118) * t141) * t140) * t113, 0, -0.2e1 * t201 + 0.2e1 * (t105 * t113 * t118 + (-t113 * t199 - t118 * t201) * t120) * t120;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end