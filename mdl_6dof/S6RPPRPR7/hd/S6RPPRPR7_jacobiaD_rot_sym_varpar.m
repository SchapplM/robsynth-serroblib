% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR7
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
%   Wie in S6RPPRPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:44
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:25
	% DurationCPUTime: 1.03s
	% Computational Cost: add. (1897->83), mult. (2191->188), div. (456->12), fcn. (2616->9), ass. (0->86)
	t104 = pkin(9) + qJ(4);
	t103 = cos(t104);
	t101 = t103 ^ 2;
	t102 = sin(t104);
	t152 = t101 / t102 ^ 2;
	t110 = cos(qJ(1));
	t128 = 0.1e1 + t152;
	t106 = t110 ^ 2;
	t96 = t106 * t152 + 0.1e1;
	t94 = 0.1e1 / t96;
	t148 = t110 * t94;
	t77 = t128 * t148;
	t167 = t110 * t77 - 0.1e1;
	t108 = cos(pkin(10));
	t109 = sin(qJ(1));
	t137 = qJD(4) * t109;
	t126 = t103 * t137;
	t107 = sin(pkin(10));
	t144 = t109 * t107;
	t145 = t108 * t110;
	t90 = t102 * t145 - t144;
	t82 = t90 * qJD(1) + t108 * t126;
	t143 = t109 * t108;
	t146 = t107 * t110;
	t88 = t102 * t143 + t146;
	t85 = 0.1e1 / t88 ^ 2;
	t166 = t82 * t85;
	t165 = t103 * t152;
	t87 = t102 * t144 - t145;
	t154 = t85 * t87;
	t83 = t87 ^ 2;
	t80 = t83 * t85 + 0.1e1;
	t78 = 0.1e1 / t80;
	t84 = 0.1e1 / t88;
	t164 = (-t107 * t84 + t108 * t154) * t78;
	t142 = t110 * t103;
	t93 = atan2(-t142, t102);
	t91 = sin(t93);
	t92 = cos(t93);
	t75 = t92 * t102 - t91 * t142;
	t72 = 0.1e1 / t75;
	t97 = 0.1e1 / t102;
	t73 = 0.1e1 / t75 ^ 2;
	t163 = t94 - 0.1e1;
	t105 = t109 ^ 2;
	t139 = qJD(1) * t110;
	t127 = t109 * t139;
	t138 = qJD(4) * t102;
	t129 = t73 * t138;
	t136 = qJD(4) * t110;
	t140 = qJD(1) * t109;
	t151 = t102 * t91;
	t68 = ((t102 * t136 + t103 * t140) * t97 + t136 * t152) * t94;
	t63 = (-t68 + t136) * t151 + (t91 * t140 + (-t110 * t68 + qJD(4)) * t92) * t103;
	t161 = t63 * t72 * t73;
	t71 = t101 * t105 * t73 + 0.1e1;
	t162 = (-t105 * t103 * t129 + (-t105 * t161 + t73 * t127) * t101) / t71 ^ 2;
	t155 = t84 * t166;
	t89 = t102 * t146 + t143;
	t81 = t89 * qJD(1) + t107 * t126;
	t160 = (t81 * t154 - t83 * t155) / t80 ^ 2;
	t69 = 0.1e1 / t71;
	t158 = t69 * t73;
	t157 = t72 * t69;
	t118 = qJD(4) * (-t103 - t165) * t97;
	t156 = (t106 * t118 - t127 * t152) / t96 ^ 2;
	t153 = t101 * t97;
	t150 = t109 * t73;
	t147 = qJD(4) * t77;
	t141 = qJD(1) * t103;
	t135 = -0.2e1 * t161;
	t134 = 0.2e1 * t160;
	t133 = t72 * t162;
	t132 = t103 * t156;
	t131 = t97 * t148;
	t130 = t103 * t163;
	t125 = t103 * t136;
	t124 = -0.2e1 * t73 * t162;
	t123 = 0.2e1 * t87 * t155;
	t122 = t101 * t131;
	t121 = t91 * t130;
	t120 = t128 * t94;
	t67 = (-t92 * t122 - t121) * t109;
	t65 = -t167 * t92 * t103 + (t110 - t77) * t151;
	t64 = -t120 * t140 + 0.2e1 * (t118 * t94 - t128 * t156) * t110;
	t1 = [t131 * t141 + (-qJD(4) * t120 - 0.2e1 * t97 * t132) * t109, 0, 0, t64, 0, 0; (t138 * t157 + (0.2e1 * t133 + (qJD(1) * t67 + t63) * t158) * t103) * t110 + (t67 * t124 * t103 + (-t67 * t129 + (t67 * t135 + ((t68 * t122 + t163 * t138 + 0.2e1 * t132) * t91 + (-t68 * t130 + (0.2e1 * t153 * t156 + (0.2e1 * t103 + t165) * t94 * qJD(4)) * t110) * t92) * t150) * t103 + (t72 + (-t110 * t121 + (t105 - t106) * t94 * t92 * t153) * t73) * t141) * t69) * t109, 0, 0, (t139 * t157 + (-0.2e1 * t133 + (-qJD(4) * t65 - t63) * t158) * t109) * t102 + (t65 * t109 * t124 + (t72 * t137 + (t109 * t135 + t73 * t139) * t65 + (((-t110 * t64 + t140 * t77) * t92 + (t167 * t68 + t136 - t147) * t91) * t103 + ((-t64 - t140) * t91 + (-t68 * t77 - qJD(4) + (t68 + t147) * t110) * t92) * t102) * t150) * t69) * t103, 0, 0; (t90 * t154 - t84 * t89) * t134 + ((-t87 * qJD(1) + t107 * t125) * t84 + t90 * t123 + (-t89 * t82 - (-t88 * qJD(1) + t108 * t125) * t87 - t90 * t81) * t85) * t78, 0, 0, t102 * t137 * t164 + (-t139 * t164 + ((-0.2e1 * t84 * t160 - t78 * t166) * t107 + (t134 * t154 + (-t81 * t85 + t123) * t78) * t108) * t109) * t103, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:25
	% DurationCPUTime: 1.02s
	% Computational Cost: add. (2429->92), mult. (2519->204), div. (480->12), fcn. (2968->9), ass. (0->95)
	t141 = pkin(9) + qJ(4);
	t137 = sin(t141);
	t132 = 0.1e1 / t137 ^ 2;
	t139 = cos(t141);
	t135 = t139 ^ 2;
	t189 = t132 * t135;
	t145 = cos(qJ(1));
	t206 = 0.2e1 * t145;
	t205 = t139 * t189;
	t182 = t145 * t139;
	t126 = atan2(-t182, t137);
	t124 = sin(t126);
	t125 = cos(t126);
	t110 = -t124 * t182 + t125 * t137;
	t107 = 0.1e1 / t110;
	t140 = pkin(10) + qJ(6);
	t136 = sin(t140);
	t138 = cos(t140);
	t144 = sin(qJ(1));
	t184 = t144 * t138;
	t121 = t136 * t145 + t137 * t184;
	t117 = 0.1e1 / t121;
	t131 = 0.1e1 / t137;
	t108 = 0.1e1 / t110 ^ 2;
	t118 = 0.1e1 / t121 ^ 2;
	t143 = t145 ^ 2;
	t129 = t143 * t189 + 0.1e1;
	t127 = 0.1e1 / t129;
	t204 = t127 - 0.1e1;
	t142 = t144 ^ 2;
	t188 = t135 * t142;
	t106 = t108 * t188 + 0.1e1;
	t179 = qJD(1) * t145;
	t161 = t135 * t144 * t179;
	t177 = qJD(4) * t139;
	t180 = qJD(1) * t144;
	t170 = t139 * t180;
	t176 = qJD(4) * t145;
	t101 = ((t137 * t176 + t170) * t131 + t176 * t189) * t127;
	t192 = t125 * t139;
	t96 = (-t101 * t145 + qJD(4)) * t192 + (t170 + (-t101 + t176) * t137) * t124;
	t202 = t107 * t108 * t96;
	t203 = 0.1e1 / t106 ^ 2 * (-t188 * t202 + (-t137 * t142 * t177 + t161) * t108);
	t164 = qJD(1) * t137 + qJD(6);
	t153 = t144 * t177 + t164 * t145;
	t165 = qJD(6) * t137 + qJD(1);
	t158 = t138 * t165;
	t102 = t153 * t136 + t144 * t158;
	t183 = t145 * t138;
	t185 = t144 * t136;
	t120 = t137 * t185 - t183;
	t116 = t120 ^ 2;
	t115 = t116 * t118 + 0.1e1;
	t194 = t118 * t120;
	t159 = t136 * t165;
	t103 = t153 * t138 - t144 * t159;
	t199 = t103 * t117 * t118;
	t201 = 0.1e1 / t115 ^ 2 * (t102 * t194 - t116 * t199);
	t200 = t101 * t139;
	t198 = t108 * t139;
	t197 = t108 * t144;
	t190 = t131 * t139;
	t156 = qJD(4) * (-t131 * t205 - t190);
	t196 = (-t132 * t161 + t143 * t156) / t129 ^ 2;
	t195 = t117 * t136;
	t193 = t120 * t138;
	t191 = t131 * t135;
	t187 = t137 * t145;
	t186 = t139 * t144;
	t168 = 0.1e1 + t189;
	t114 = t168 * t145 * t127;
	t181 = -t114 + t145;
	t178 = qJD(4) * t137;
	t175 = -0.2e1 * t202;
	t174 = 0.2e1 * t201;
	t173 = t139 * t203;
	t172 = t139 * t196;
	t171 = t127 * t191;
	t169 = t114 * t145 - 0.1e1;
	t167 = 0.2e1 * t120 * t199;
	t166 = t196 * t206;
	t163 = t145 * t171;
	t162 = t204 * t139 * t124;
	t160 = t168 * t144;
	t157 = t118 * t193 - t195;
	t155 = t157 * t144;
	t154 = t139 * t176 - t164 * t144;
	t123 = t137 * t183 - t185;
	t122 = t136 * t187 + t184;
	t112 = 0.1e1 / t115;
	t104 = 0.1e1 / t106;
	t100 = (-t125 * t163 - t162) * t144;
	t99 = t124 * t187 + t192 + (-t124 * t137 - t125 * t182) * t114;
	t97 = -t168 * t166 + (-qJD(1) * t160 + t156 * t206) * t127;
	t1 = [-0.2e1 * t144 * t131 * t172 + (-qJD(4) * t160 + t179 * t190) * t127, 0, 0, t97, 0, 0; (0.2e1 * t107 * t173 + (t107 * t178 + (qJD(1) * t100 + t96) * t198) * t104) * t145 + (-0.2e1 * t108 * t173 * t100 + (((t101 * t163 + t204 * t178 + 0.2e1 * t172) * t124 + (t166 * t191 + t200 + (-t200 + (0.2e1 * t139 + t205) * t176) * t127) * t125) * t108 * t186 + (-t108 * t178 + t139 * t175) * t100 + (t107 + ((t142 - t143) * t125 * t171 - t145 * t162) * t108) * t139 * qJD(1)) * t104) * t144, 0, 0, 0.2e1 * (-t107 * t137 - t99 * t198) * t144 * t203 + ((t107 * t179 + (-qJD(4) * t99 - t96) * t197) * t137 + ((qJD(4) * t107 + t99 * t175) * t144 + (t99 * t179 + ((t114 * t180 - t145 * t97) * t125 + (t181 * qJD(4) + t169 * t101) * t124) * t186) * t108 + ((-t97 - t180) * t124 + (t169 * qJD(4) + t181 * t101) * t125) * t137 * t197) * t139) * t104, 0, 0; (-t117 * t122 + t123 * t194) * t174 + (t123 * t167 + t145 * t117 * t158 + t154 * t195 + (t145 * t120 * t159 - t123 * t102 - t122 * t103 - t154 * t193) * t118) * t112, 0, 0, t139 * t155 * t174 + (t155 * t178 + (-t157 * t179 + ((qJD(6) * t117 + t167) * t138 + (-t102 * t138 + (qJD(6) * t120 - t103) * t136) * t118) * t144) * t139) * t112, 0, -0.2e1 * t201 + 0.2e1 * (t102 * t112 * t118 + (-t112 * t199 - t118 * t201) * t120) * t120;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end