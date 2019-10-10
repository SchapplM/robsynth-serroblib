% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR2
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
%   Wie in S6RPRPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:17
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:01
	% DurationCPUTime: 0.80s
	% Computational Cost: add. (2561->72), mult. (1839->154), div. (436->14), fcn. (2165->7), ass. (0->77)
	t102 = qJ(1) + pkin(9);
	t100 = cos(t102);
	t129 = qJD(1) * t100;
	t98 = sin(t102);
	t153 = 0.2e1 * t98;
	t87 = t98 ^ 2;
	t88 = 0.1e1 / t98;
	t96 = t100 ^ 2;
	t151 = (0.1e1 + 0.1e1 / t87 * t96) * t88 * t129;
	t101 = qJ(3) + pkin(10);
	t97 = sin(t101);
	t136 = t98 * t97;
	t99 = cos(t101);
	t78 = atan2(-t136, -t99);
	t76 = sin(t78);
	t124 = t76 * t136;
	t77 = cos(t78);
	t72 = -t77 * t99 - t124;
	t69 = 0.1e1 / t72;
	t92 = 0.1e1 / t99;
	t150 = -0.2e1 * t97;
	t70 = 0.1e1 / t72 ^ 2;
	t86 = t97 ^ 2;
	t93 = 0.1e1 / t99 ^ 2;
	t139 = t86 * t93;
	t83 = t87 * t139 + 0.1e1;
	t79 = 0.1e1 / t83;
	t149 = t79 - 0.1e1;
	t119 = t97 * t129;
	t131 = qJD(3) * t98;
	t141 = t77 * t97;
	t130 = qJD(3) * t99;
	t65 = (-(-t98 * t130 - t119) * t92 + t131 * t139) * t79;
	t61 = (-t65 * t98 + qJD(3)) * t141 + (-t119 + (t65 - t131) * t99) * t76;
	t148 = t61 * t69 * t70;
	t147 = t65 * t76;
	t146 = t65 * t97;
	t145 = t70 * t97;
	t85 = t97 * t86;
	t91 = t99 ^ 2;
	t109 = qJD(3) * (t85 / t91 + t97) * t92;
	t113 = t86 * t98 * t129;
	t144 = (t87 * t109 + t93 * t113) / t83 ^ 2;
	t121 = 0.1e1 + t139;
	t112 = t121 * t79;
	t75 = t98 * t112;
	t143 = t75 * t98;
	t142 = t76 * t99;
	t140 = t86 * t92;
	t138 = t86 * t96;
	t89 = 0.1e1 / t98 ^ 2;
	t137 = t89 * t96;
	t135 = t100 * t70;
	t134 = t100 * t97;
	t133 = qJD(1) * t97;
	t132 = qJD(1) * t98;
	t128 = qJD(3) * t100;
	t114 = t96 * t97 * t130;
	t68 = t70 * t138 + 0.1e1;
	t127 = 0.2e1 * (-t138 * t148 + (-t113 + t114) * t70) / t68 ^ 2;
	t126 = 0.2e1 * t148;
	t84 = t91 * t137 + 0.1e1;
	t125 = 0.2e1 * (-t89 * t114 - t91 * t151) / t84 ^ 2;
	t123 = t79 * t92 * t98;
	t122 = t70 * t134;
	t120 = 0.1e1 + t137;
	t118 = t97 * t127;
	t117 = t144 * t150;
	t116 = t144 * t153;
	t115 = t86 * t123;
	t111 = t120 * t97;
	t81 = 0.1e1 / t84;
	t66 = 0.1e1 / t68;
	t64 = (t149 * t97 * t76 - t77 * t115) * t100;
	t63 = -t98 * t142 + t141 + (-t77 * t136 + t142) * t75;
	t62 = -t121 * t116 + (t109 * t153 + t121 * t129) * t79;
	t1 = [-t123 * t133 + (qJD(3) * t112 + t92 * t117) * t100, 0, t62, 0, 0, 0; (t69 * t118 + (-t69 * t130 + (qJD(1) * t64 + t61) * t145) * t66) * t98 + (t70 * t118 * t64 + (-((t65 * t115 + t149 * t130 + t117) * t76 + (t116 * t140 - t146 + (t146 + (-t85 * t93 + t150) * t131) * t79) * t77) * t122 + (t97 * t126 - t70 * t130) * t64 + (-t69 + ((-t87 + t96) * t79 * t77 * t140 + t149 * t124) * t70) * t133) * t66) * t100, 0, (t63 * t145 - t69 * t99) * t100 * t127 + ((-t69 * t132 + (-qJD(3) * t63 - t61) * t135) * t99 + (-t69 * t128 - (-t62 * t77 * t98 + t76 * t131 + t143 * t147 - t147 + (-qJD(3) * t76 - t129 * t77) * t75) * t122 + (t100 * t126 + t70 * t132) * t63 - ((t62 - t129) * t76 + ((0.1e1 - t143) * qJD(3) + (t75 - t98) * t65) * t77) * t99 * t135) * t97) * t66, 0, 0, 0; t120 * t99 * t125 + (qJD(3) * t111 + 0.2e1 * t99 * t151) * t81, 0, t88 * t125 * t134 + (-t88 * t99 * t128 + qJD(1) * t111) * t81, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:01
	% DurationCPUTime: 1.06s
	% Computational Cost: add. (3055->93), mult. (2519->205), div. (480->12), fcn. (2968->9), ass. (0->93)
	t145 = qJ(3) + pkin(10);
	t141 = sin(t145);
	t135 = 0.1e1 / t141 ^ 2;
	t143 = cos(t145);
	t139 = t143 ^ 2;
	t193 = t135 * t139;
	t209 = t143 * t193;
	t146 = qJ(1) + pkin(9);
	t142 = sin(t146);
	t168 = 0.1e1 + t193;
	t208 = t142 * t168;
	t137 = t142 ^ 2;
	t132 = t137 * t193 + 0.1e1;
	t130 = 0.1e1 / t132;
	t134 = 0.1e1 / t141;
	t144 = cos(t146);
	t183 = qJD(1) * t144;
	t171 = t143 * t183;
	t181 = qJD(3) * t142;
	t104 = ((t141 * t181 - t171) * t134 + t181 * t193) * t130;
	t207 = -t104 + t181;
	t147 = sin(qJ(6));
	t148 = cos(qJ(6));
	t165 = qJD(6) * t141 + qJD(1);
	t180 = qJD(3) * t143;
	t206 = t165 * t147 - t148 * t180;
	t205 = t147 * t180 + t165 * t148;
	t190 = t142 * t143;
	t129 = atan2(-t190, t141);
	t128 = cos(t129);
	t127 = sin(t129);
	t173 = t127 * t190;
	t115 = t128 * t141 - t173;
	t111 = 0.1e1 / t115;
	t187 = t144 * t147;
	t188 = t142 * t148;
	t124 = t141 * t187 + t188;
	t120 = 0.1e1 / t124;
	t112 = 0.1e1 / t115 ^ 2;
	t121 = 0.1e1 / t124 ^ 2;
	t204 = t130 - 0.1e1;
	t194 = t128 * t143;
	t99 = (-t104 * t142 + qJD(3)) * t194 + (t207 * t141 - t171) * t127;
	t203 = t111 * t112 * t99;
	t164 = qJD(1) * t141 + qJD(6);
	t159 = t164 * t148;
	t108 = t142 * t159 + t206 * t144;
	t186 = t144 * t148;
	t189 = t142 * t147;
	t123 = -t141 * t186 + t189;
	t119 = t123 ^ 2;
	t118 = t119 * t121 + 0.1e1;
	t196 = t121 * t123;
	t160 = t164 * t147;
	t109 = -t142 * t160 + t205 * t144;
	t200 = t109 * t120 * t121;
	t202 = (t108 * t196 - t119 * t200) / t118 ^ 2;
	t201 = t104 * t143;
	t157 = qJD(3) * (-t143 - t209) * t134;
	t191 = t139 * t142;
	t162 = t183 * t191;
	t199 = (t135 * t162 + t137 * t157) / t132 ^ 2;
	t198 = t112 * t143;
	t197 = t112 * t144;
	t195 = t127 * t142;
	t140 = t144 ^ 2;
	t192 = t139 * t140;
	t185 = qJD(1) * t142;
	t184 = qJD(1) * t143;
	t182 = qJD(3) * t141;
	t107 = t112 * t192 + 0.1e1;
	t179 = 0.2e1 / t107 ^ 2 * (-t192 * t203 + (-t140 * t141 * t180 - t162) * t112);
	t178 = 0.2e1 * t203;
	t177 = 0.2e1 * t202;
	t176 = -0.2e1 * t199;
	t175 = t143 * t199;
	t174 = t143 * t197;
	t172 = t134 * t191;
	t167 = t143 * t179;
	t166 = 0.2e1 * t123 * t200;
	t163 = t130 * t172;
	t161 = t168 * t144;
	t158 = t120 * t148 + t147 * t196;
	t156 = t158 * t144;
	t126 = -t141 * t189 + t186;
	t125 = t141 * t188 + t187;
	t116 = 0.1e1 / t118;
	t114 = t130 * t208;
	t105 = 0.1e1 / t107;
	t103 = (t204 * t143 * t127 + t128 * t163) * t144;
	t101 = t141 * t195 + t194 + (-t127 * t141 - t128 * t190) * t114;
	t100 = t176 * t208 + (qJD(1) * t161 + 0.2e1 * t142 * t157) * t130;
	t1 = [0.2e1 * t134 * t144 * t175 + (t134 * t142 * t184 + qJD(3) * t161) * t130, 0, t100, 0, 0, 0; (t111 * t167 + (t111 * t182 + (qJD(1) * t103 + t99) * t198) * t105) * t142 + (t112 * t167 * t103 + (-((-t104 * t163 - t204 * t182 - 0.2e1 * t175) * t127 + (t172 * t176 - t201 + (t201 + (-0.2e1 * t143 - t209) * t181) * t130) * t128) * t174 + (t112 * t182 + t143 * t178) * t103 + (-t111 + ((t137 - t140) * t139 * t134 * t130 * t128 + t204 * t173) * t112) * t184) * t105) * t144, 0, (t101 * t198 + t111 * t141) * t144 * t179 + ((t111 * t185 + (qJD(3) * t101 + t99) * t197) * t141 + (-t144 * qJD(3) * t111 - (-t100 * t128 * t142 + t207 * t127 + (-qJD(3) * t127 + t104 * t195 - t128 * t183) * t114) * t174 + (t112 * t185 + t144 * t178) * t101 - ((-t100 + t183) * t127 + ((t114 * t142 - 0.1e1) * qJD(3) + (-t114 + t142) * t104) * t128) * t141 * t197) * t143) * t105, 0, 0, 0; (-t120 * t125 + t126 * t196) * t177 + (t126 * t166 + (-t126 * t108 - t125 * t109 + (t205 * t142 + t144 * t160) * t123) * t121 + (-t206 * t142 + t144 * t159) * t120) * t116, 0, t143 * t156 * t177 + (t156 * t182 + (t158 * t185 + ((qJD(6) * t120 + t166) * t147 + (-t108 * t147 + (-qJD(6) * t123 + t109) * t148) * t121) * t144) * t143) * t116, 0, 0, -0.2e1 * t202 + 0.2e1 * (t108 * t116 * t121 + (-t116 * t200 - t121 * t202) * t123) * t123;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end