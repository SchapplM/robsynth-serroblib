% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR6
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
%   Wie in S6RPRPPR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:24
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:01
	% EndTime: 2019-10-10 00:24:02
	% DurationCPUTime: 0.97s
	% Computational Cost: add. (1897->83), mult. (2191->188), div. (456->12), fcn. (2616->9), ass. (0->86)
	t107 = qJ(3) + pkin(9);
	t106 = cos(t107);
	t104 = t106 ^ 2;
	t105 = sin(t107);
	t150 = 0.1e1 / t105 ^ 2 * t104;
	t113 = cos(qJ(1));
	t126 = 0.1e1 + t150;
	t109 = t113 ^ 2;
	t99 = t109 * t150 + 0.1e1;
	t97 = 0.1e1 / t99;
	t153 = t113 * t97;
	t80 = t126 * t153;
	t170 = t113 * t80 - 0.1e1;
	t111 = cos(pkin(10));
	t112 = sin(qJ(1));
	t140 = qJD(3) * t112;
	t130 = t106 * t140;
	t110 = sin(pkin(10));
	t147 = t112 * t110;
	t148 = t111 * t113;
	t93 = t105 * t148 - t147;
	t85 = t93 * qJD(1) + t111 * t130;
	t146 = t112 * t111;
	t149 = t110 * t113;
	t91 = t105 * t146 + t149;
	t88 = 0.1e1 / t91 ^ 2;
	t169 = t85 * t88;
	t90 = t105 * t147 - t148;
	t157 = t88 * t90;
	t86 = t90 ^ 2;
	t83 = t86 * t88 + 0.1e1;
	t81 = 0.1e1 / t83;
	t87 = 0.1e1 / t91;
	t168 = (-t110 * t87 + t111 * t157) * t81;
	t167 = t106 * t150;
	t145 = t113 * t106;
	t96 = atan2(-t145, t105);
	t94 = sin(t96);
	t95 = cos(t96);
	t78 = t95 * t105 - t94 * t145;
	t75 = 0.1e1 / t78;
	t100 = 0.1e1 / t105;
	t76 = 0.1e1 / t78 ^ 2;
	t166 = t97 - 0.1e1;
	t108 = t112 ^ 2;
	t142 = qJD(1) * t113;
	t131 = t112 * t142;
	t141 = qJD(3) * t105;
	t132 = t76 * t141;
	t139 = qJD(3) * t113;
	t143 = qJD(1) * t112;
	t156 = t105 * t94;
	t71 = ((t105 * t139 + t106 * t143) * t100 + t139 * t150) * t97;
	t66 = (-t71 + t139) * t156 + (t94 * t143 + (-t113 * t71 + qJD(3)) * t95) * t106;
	t164 = t66 * t75 * t76;
	t74 = t104 * t108 * t76 + 0.1e1;
	t165 = (-t108 * t106 * t132 + (-t108 * t164 + t76 * t131) * t104) / t74 ^ 2;
	t158 = t87 * t169;
	t92 = t105 * t149 + t146;
	t84 = t92 * qJD(1) + t110 * t130;
	t163 = (t84 * t157 - t86 * t158) / t83 ^ 2;
	t72 = 0.1e1 / t74;
	t161 = t72 * t76;
	t160 = t75 * t72;
	t121 = qJD(3) * (-t106 - t167) * t100;
	t159 = (t109 * t121 - t131 * t150) / t99 ^ 2;
	t155 = t112 * t76;
	t152 = qJD(3) * t80;
	t151 = t100 * t104;
	t144 = qJD(1) * t106;
	t138 = -0.2e1 * t164;
	t137 = 0.2e1 * t163;
	t136 = t75 * t165;
	t135 = t106 * t159;
	t134 = t100 * t153;
	t133 = t106 * t166;
	t129 = t106 * t139;
	t128 = -0.2e1 * t76 * t165;
	t127 = 0.2e1 * t90 * t158;
	t125 = t94 * t133;
	t124 = t104 * t134;
	t123 = t126 * t97;
	t70 = (-t95 * t124 - t125) * t112;
	t68 = -t170 * t95 * t106 + (t113 - t80) * t156;
	t67 = -t123 * t143 + 0.2e1 * (t121 * t97 - t126 * t159) * t113;
	t1 = [t134 * t144 + (-qJD(3) * t123 - 0.2e1 * t100 * t135) * t112, 0, t67, 0, 0, 0; (t141 * t160 + (0.2e1 * t136 + (qJD(1) * t70 + t66) * t161) * t106) * t113 + (t70 * t128 * t106 + (-t70 * t132 + (t70 * t138 + ((t71 * t124 + t166 * t141 + 0.2e1 * t135) * t94 + (-t71 * t133 + (0.2e1 * t151 * t159 + (0.2e1 * t106 + t167) * t97 * qJD(3)) * t113) * t95) * t155) * t106 + (t75 + (-t113 * t125 + (t108 - t109) * t97 * t95 * t151) * t76) * t144) * t72) * t112, 0, (t142 * t160 + (-0.2e1 * t136 + (-qJD(3) * t68 - t66) * t161) * t112) * t105 + (t68 * t112 * t128 + (t75 * t140 + (t112 * t138 + t76 * t142) * t68 + (((-t113 * t67 + t143 * t80) * t95 + (t170 * t71 + t139 - t152) * t94) * t106 + ((-t67 - t143) * t94 + (-t71 * t80 - qJD(3) + (t71 + t152) * t113) * t95) * t105) * t155) * t72) * t106, 0, 0, 0; (t93 * t157 - t87 * t92) * t137 + ((-t90 * qJD(1) + t110 * t129) * t87 + t93 * t127 + (-t92 * t85 - (-t91 * qJD(1) + t111 * t129) * t90 - t93 * t84) * t88) * t81, 0, t105 * t140 * t168 + (-t142 * t168 + ((-0.2e1 * t87 * t163 - t81 * t169) * t110 + (t137 * t157 + (-t84 * t88 + t127) * t81) * t111) * t112) * t106, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:01
	% EndTime: 2019-10-10 00:24:02
	% DurationCPUTime: 1.01s
	% Computational Cost: add. (2429->94), mult. (2519->204), div. (480->12), fcn. (2968->9), ass. (0->95)
	t144 = qJ(3) + pkin(9);
	t140 = sin(t144);
	t135 = 0.1e1 / t140 ^ 2;
	t142 = cos(t144);
	t138 = t142 ^ 2;
	t190 = t135 * t138;
	t148 = cos(qJ(1));
	t209 = 0.2e1 * t148;
	t208 = t142 * t190;
	t184 = t148 * t142;
	t129 = atan2(-t184, t140);
	t127 = sin(t129);
	t128 = cos(t129);
	t113 = -t127 * t184 + t128 * t140;
	t110 = 0.1e1 / t113;
	t143 = pkin(10) + qJ(6);
	t139 = sin(t143);
	t141 = cos(t143);
	t147 = sin(qJ(1));
	t186 = t147 * t141;
	t124 = t139 * t148 + t140 * t186;
	t120 = 0.1e1 / t124;
	t134 = 0.1e1 / t140;
	t111 = 0.1e1 / t113 ^ 2;
	t121 = 0.1e1 / t124 ^ 2;
	t146 = t148 ^ 2;
	t132 = t146 * t190 + 0.1e1;
	t130 = 0.1e1 / t132;
	t207 = t130 - 0.1e1;
	t145 = t147 ^ 2;
	t189 = t138 * t145;
	t109 = t111 * t189 + 0.1e1;
	t182 = qJD(1) * t148;
	t164 = t138 * t147 * t182;
	t181 = qJD(3) * t140;
	t183 = qJD(1) * t147;
	t172 = t142 * t183;
	t179 = qJD(3) * t148;
	t104 = ((t140 * t179 + t172) * t134 + t179 * t190) * t130;
	t193 = t128 * t142;
	t99 = (-t104 * t148 + qJD(3)) * t193 + (t172 + (-t104 + t179) * t140) * t127;
	t205 = t110 * t111 * t99;
	t206 = 0.1e1 / t109 ^ 2 * (-t189 * t205 + (-t142 * t145 * t181 + t164) * t111);
	t167 = qJD(1) * t140 + qJD(6);
	t180 = qJD(3) * t147;
	t156 = t142 * t180 + t167 * t148;
	t168 = qJD(6) * t140 + qJD(1);
	t161 = t141 * t168;
	t105 = t156 * t139 + t147 * t161;
	t185 = t148 * t141;
	t187 = t147 * t139;
	t123 = t140 * t187 - t185;
	t119 = t123 ^ 2;
	t118 = t119 * t121 + 0.1e1;
	t195 = t121 * t123;
	t162 = t139 * t168;
	t106 = t156 * t141 - t147 * t162;
	t201 = t106 * t120 * t121;
	t204 = (t105 * t195 - t119 * t201) / t118 ^ 2;
	t203 = t104 * t127;
	t202 = t104 * t142;
	t200 = t111 * t142;
	t199 = t111 * t147;
	t191 = t134 * t142;
	t159 = qJD(3) * (-t134 * t208 - t191);
	t198 = (-t135 * t164 + t146 * t159) / t132 ^ 2;
	t171 = 0.1e1 + t190;
	t117 = t171 * t148 * t130;
	t197 = t117 * t148;
	t196 = t120 * t139;
	t194 = t123 * t141;
	t192 = t134 * t138;
	t188 = t140 * t148;
	t178 = -0.2e1 * t205;
	t177 = 0.2e1 * t204;
	t176 = t142 * t206;
	t175 = t142 * t199;
	t174 = t142 * t198;
	t173 = t130 * t192;
	t170 = 0.2e1 * t123 * t201;
	t169 = t198 * t209;
	t166 = t148 * t173;
	t165 = t207 * t142 * t127;
	t163 = t171 * t147;
	t160 = t121 * t194 - t196;
	t158 = t160 * t147;
	t157 = t142 * t179 - t167 * t147;
	t126 = t140 * t185 - t187;
	t125 = t139 * t188 + t186;
	t115 = 0.1e1 / t118;
	t107 = 0.1e1 / t109;
	t103 = (-t128 * t166 - t165) * t147;
	t102 = t127 * t188 + t193 + (-t127 * t140 - t128 * t184) * t117;
	t100 = -t171 * t169 + (-qJD(1) * t163 + t159 * t209) * t130;
	t1 = [-0.2e1 * t147 * t134 * t174 + (-qJD(3) * t163 + t182 * t191) * t130, 0, t100, 0, 0, 0; (0.2e1 * t110 * t176 + (t110 * t181 + (qJD(1) * t103 + t99) * t200) * t107) * t148 + (-0.2e1 * t111 * t176 * t103 + (((t104 * t166 + t207 * t181 + 0.2e1 * t174) * t127 + (t169 * t192 + t202 + (-t202 + (0.2e1 * t142 + t208) * t179) * t130) * t128) * t175 + (-t111 * t181 + t142 * t178) * t103 + (t110 + ((t145 - t146) * t128 * t173 - t148 * t165) * t111) * t142 * qJD(1)) * t107) * t147, 0, 0.2e1 * (-t102 * t200 - t110 * t140) * t147 * t206 + ((t110 * t182 + (-qJD(3) * t102 - t99) * t199) * t140 + (t110 * t180 + (-t100 * t128 * t148 + t127 * t179 + t197 * t203 - t203 + (-qJD(3) * t127 + t128 * t183) * t117) * t175 + (t111 * t182 + t147 * t178) * t102 + ((-t100 - t183) * t127 + ((-0.1e1 + t197) * qJD(3) + (-t117 + t148) * t104) * t128) * t140 * t199) * t142) * t107, 0, 0, 0; (-t120 * t125 + t126 * t195) * t177 + (t126 * t170 + t148 * t120 * t161 + t157 * t196 + (t148 * t123 * t162 - t126 * t105 - t125 * t106 - t157 * t194) * t121) * t115, 0, t142 * t158 * t177 + (t158 * t181 + (-t160 * t182 + ((qJD(6) * t120 + t170) * t141 + (-t105 * t141 + (qJD(6) * t123 - t106) * t139) * t121) * t147) * t142) * t115, 0, 0, -0.2e1 * t204 + 0.2e1 * (t105 * t115 * t121 + (-t115 * t201 - t121 * t204) * t123) * t123;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end