% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPPR6
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
%   Wie in S5RRPPR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRPPR6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:54
	% EndTime: 2019-12-29 18:19:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:49
	% EndTime: 2019-12-29 18:19:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:49
	% EndTime: 2019-12-29 18:19:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:49
	% EndTime: 2019-12-29 18:19:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:54
	% EndTime: 2019-12-29 18:19:55
	% DurationCPUTime: 1.49s
	% Computational Cost: add. (2086->84), mult. (2191->187), div. (456->12), fcn. (2616->9), ass. (0->87)
	t107 = qJ(2) + pkin(8);
	t105 = sin(t107);
	t101 = t105 ^ 2;
	t106 = cos(t107);
	t150 = t101 / t106 ^ 2;
	t112 = sin(qJ(1));
	t126 = 0.1e1 + t150;
	t108 = t112 ^ 2;
	t99 = t108 * t150 + 0.1e1;
	t97 = 0.1e1 / t99;
	t123 = t126 * t97;
	t80 = t112 * t123;
	t171 = t112 * t80 - 0.1e1;
	t111 = cos(pkin(9));
	t113 = cos(qJ(1));
	t139 = qJD(2) * t113;
	t128 = t105 * t139;
	t145 = t112 * t111;
	t110 = sin(pkin(9));
	t149 = t110 * t113;
	t91 = -t106 * t145 + t149;
	t85 = t91 * qJD(1) - t111 * t128;
	t146 = t112 * t110;
	t148 = t111 * t113;
	t93 = t106 * t148 + t146;
	t88 = 0.1e1 / t93 ^ 2;
	t170 = t85 * t88;
	t92 = t106 * t149 - t145;
	t157 = t88 * t92;
	t86 = t92 ^ 2;
	t83 = t86 * t88 + 0.1e1;
	t81 = 0.1e1 / t83;
	t87 = 0.1e1 / t93;
	t169 = (-t110 * t87 + t111 * t157) * t81;
	t168 = t105 * t150;
	t147 = t112 * t105;
	t96 = atan2(-t147, -t106);
	t94 = sin(t96);
	t133 = t94 * t147;
	t95 = cos(t96);
	t78 = -t106 * t95 - t133;
	t75 = 0.1e1 / t78;
	t102 = 0.1e1 / t106;
	t76 = 0.1e1 / t78 ^ 2;
	t167 = 0.2e1 * t105;
	t166 = t97 - 0.1e1;
	t109 = t113 ^ 2;
	t142 = qJD(1) * t113;
	t130 = t112 * t142;
	t141 = qJD(2) * t106;
	t131 = t76 * t141;
	t140 = qJD(2) * t112;
	t154 = t106 * t94;
	t71 = (-(-t105 * t142 - t106 * t140) * t102 + t140 * t150) * t97;
	t66 = (t71 - t140) * t154 + (-t94 * t142 + (-t112 * t71 + qJD(2)) * t95) * t105;
	t164 = t66 * t75 * t76;
	t156 = t101 * t76;
	t74 = t109 * t156 + 0.1e1;
	t165 = (t109 * t105 * t131 + (-t109 * t164 - t76 * t130) * t101) / t74 ^ 2;
	t158 = t87 * t170;
	t90 = -t106 * t146 - t148;
	t84 = t90 * qJD(1) - t110 * t128;
	t163 = (t84 * t157 - t86 * t158) / t83 ^ 2;
	t72 = 0.1e1 / t74;
	t161 = t72 * t76;
	t160 = t75 * t72;
	t121 = qJD(2) * (t105 + t168) * t102;
	t159 = (t108 * t121 + t130 * t150) / t99 ^ 2;
	t155 = t102 * t97;
	t152 = t113 * t76;
	t151 = qJD(2) * t80;
	t144 = qJD(1) * t105;
	t143 = qJD(1) * t112;
	t138 = 0.2e1 * t164;
	t137 = 0.2e1 * t163;
	t136 = t75 * t165;
	t135 = t92 * t158;
	t134 = t112 * t155;
	t132 = t166 * t105;
	t129 = t105 * t140;
	t127 = 0.2e1 * t76 * t165;
	t125 = -0.2e1 * t102 * t159;
	t124 = t101 * t134;
	t70 = (-t95 * t124 + t94 * t132) * t113;
	t68 = (-t112 + t80) * t154 - t171 * t95 * t105;
	t67 = t123 * t142 + 0.2e1 * (t121 * t97 - t126 * t159) * t112;
	t1 = [-t134 * t144 + (qJD(2) * t123 + t105 * t125) * t113, t67, 0, 0, 0; (-t141 * t160 + (0.2e1 * t136 + (qJD(1) * t70 + t66) * t161) * t105) * t112 + (t70 * t127 * t105 + (-t70 * t131 + (t70 * t138 + ((-t71 * t124 - t166 * t141 + t159 * t167) * t94 + (-t71 * t132 + (t101 * t125 + (t167 + t168) * t97 * qJD(2)) * t112) * t95) * t152) * t105 + (-t75 + t166 * t76 * t133 - (t108 - t109) * t95 * t155 * t156) * t144) * t72) * t113, (-t143 * t160 + (-0.2e1 * t136 + (-qJD(2) * t68 - t66) * t161) * t113) * t106 + (t68 * t113 * t127 + (-t75 * t139 + (t113 * t138 + t76 * t143) * t68 + (-((-t112 * t67 - t142 * t80) * t95 + (t171 * t71 + t140 - t151) * t94) * t105 - ((t67 - t142) * t94 + (t71 * t80 + qJD(2) + (-t71 - t151) * t112) * t95) * t106) * t152) * t72) * t105, 0, 0, 0; (t91 * t157 - t87 * t90) * t137 + ((-t92 * qJD(1) + t110 * t129) * t87 + 0.2e1 * t91 * t135 + (-t90 * t85 - (-t93 * qJD(1) + t111 * t129) * t92 - t91 * t84) * t88) * t81, t106 * t139 * t169 + (-t143 * t169 + ((t87 * t137 + t81 * t170) * t110 + (-0.2e1 * t157 * t163 + (t84 * t88 - 0.2e1 * t135) * t81) * t111) * t113) * t105, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:49
	% EndTime: 2019-12-29 18:19:51
	% DurationCPUTime: 1.66s
	% Computational Cost: add. (2618->95), mult. (2519->203), div. (480->12), fcn. (2968->9), ass. (0->93)
	t145 = qJ(2) + pkin(8);
	t141 = sin(t145);
	t136 = t141 ^ 2;
	t143 = cos(t145);
	t138 = 0.1e1 / t143 ^ 2;
	t193 = t136 * t138;
	t148 = sin(qJ(1));
	t210 = 0.2e1 * t148;
	t209 = t141 * t193;
	t146 = t148 ^ 2;
	t131 = t146 * t193 + 0.1e1;
	t129 = 0.1e1 / t131;
	t137 = 0.1e1 / t143;
	t149 = cos(qJ(1));
	t183 = qJD(1) * t149;
	t171 = t141 * t183;
	t181 = qJD(2) * t148;
	t103 = (-(-t143 * t181 - t171) * t137 + t181 * t193) * t129;
	t208 = t103 - t181;
	t144 = pkin(9) + qJ(5);
	t142 = cos(t144);
	t140 = sin(t144);
	t187 = t148 * t140;
	t188 = t143 * t149;
	t125 = t142 * t188 + t187;
	t186 = t148 * t141;
	t128 = atan2(-t186, -t143);
	t127 = cos(t128);
	t126 = sin(t128);
	t174 = t126 * t186;
	t112 = -t127 * t143 - t174;
	t109 = 0.1e1 / t112;
	t119 = 0.1e1 / t125;
	t110 = 0.1e1 / t112 ^ 2;
	t120 = 0.1e1 / t125 ^ 2;
	t207 = -0.2e1 * t141;
	t206 = t129 - 0.1e1;
	t195 = t127 * t141;
	t98 = (-t103 * t148 + qJD(2)) * t195 + (t143 * t208 - t171) * t126;
	t205 = t109 * t110 * t98;
	t159 = t142 * t149 + t143 * t187;
	t180 = qJD(2) * t149;
	t170 = t141 * t180;
	t104 = t159 * qJD(1) - t125 * qJD(5) + t140 * t170;
	t185 = t148 * t142;
	t124 = t140 * t188 - t185;
	t118 = t124 ^ 2;
	t117 = t118 * t120 + 0.1e1;
	t198 = t120 * t124;
	t164 = -qJD(1) * t143 + qJD(5);
	t165 = qJD(5) * t143 - qJD(1);
	t190 = t140 * t149;
	t105 = -t165 * t190 + (t164 * t148 - t170) * t142;
	t202 = t105 * t119 * t120;
	t204 = (-t104 * t198 - t118 * t202) / t117 ^ 2;
	t203 = t103 * t141;
	t201 = t110 * t141;
	t191 = t137 * t141;
	t158 = qJD(2) * (t137 * t209 + t191);
	t162 = t136 * t148 * t183;
	t200 = (t138 * t162 + t146 * t158) / t131 ^ 2;
	t199 = t119 * t140;
	t197 = t124 * t142;
	t196 = t126 * t148;
	t194 = t136 * t137;
	t147 = t149 ^ 2;
	t192 = t136 * t147;
	t189 = t141 * t149;
	t184 = qJD(1) * t148;
	t182 = qJD(2) * t143;
	t108 = t110 * t192 + 0.1e1;
	t179 = 0.2e1 / t108 ^ 2 * (-t192 * t205 + (t141 * t147 * t182 - t162) * t110);
	t178 = 0.2e1 * t205;
	t177 = -0.2e1 * t204;
	t176 = t124 * t202;
	t175 = t110 * t189;
	t173 = t129 * t194;
	t169 = 0.1e1 + t193;
	t168 = t141 * t179;
	t167 = t200 * t207;
	t166 = t200 * t210;
	t163 = t148 * t173;
	t161 = t169 * t149;
	t160 = t120 * t197 - t199;
	t157 = t141 * t181 + t164 * t149;
	t123 = -t143 * t185 + t190;
	t116 = t169 * t148 * t129;
	t114 = 0.1e1 / t117;
	t106 = 0.1e1 / t108;
	t102 = (t206 * t141 * t126 - t127 * t163) * t149;
	t101 = -t143 * t196 + t195 + (t126 * t143 - t127 * t186) * t116;
	t99 = -t169 * t166 + (qJD(1) * t161 + t158 * t210) * t129;
	t1 = [t137 * t149 * t167 + (qJD(2) * t161 - t184 * t191) * t129, t99, 0, 0, 0; (t109 * t168 + (-t109 * t182 + (qJD(1) * t102 + t98) * t201) * t106) * t148 + (t110 * t168 * t102 + (-((t103 * t163 + t206 * t182 + t167) * t126 + (t166 * t194 - t203 + (t203 + (t207 - t209) * t181) * t129) * t127) * t175 + (-t110 * t182 + t141 * t178) * t102 + (-t109 + ((-t146 + t147) * t127 * t173 + t206 * t174) * t110) * t141 * qJD(1)) * t106) * t149, (t101 * t201 - t109 * t143) * t149 * t179 + ((-t109 * t184 + (-qJD(2) * t101 - t98) * t149 * t110) * t143 + (-t109 * t180 - (-t127 * t148 * t99 - t208 * t126 + (-qJD(2) * t126 + t103 * t196 - t127 * t183) * t116) * t175 + (t110 * t184 + t149 * t178) * t101 - ((t99 - t183) * t126 + ((-t116 * t148 + 0.1e1) * qJD(2) + (t116 - t148) * t103) * t127) * t110 * t188) * t141) * t106, 0, 0, 0; 0.2e1 * (t119 * t159 + t123 * t198) * t204 + (0.2e1 * t123 * t176 - t165 * t119 * t185 + t157 * t199 + (-t165 * t124 * t187 + t123 * t104 + t105 * t159 - t157 * t197) * t120) * t114, t160 * t177 * t189 + (t160 * t143 * t180 + (-t160 * t184 + ((-qJD(5) * t119 - 0.2e1 * t176) * t142 + (-t104 * t142 + (-qJD(5) * t124 + t105) * t140) * t120) * t149) * t141) * t114, 0, 0, t177 + 0.2e1 * (-t104 * t114 * t120 + (-t114 * t202 - t120 * t204) * t124) * t124;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end