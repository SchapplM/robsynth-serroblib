% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR1
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
%   Wie in S6RPRPPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:15
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:14
	% EndTime: 2019-10-10 00:15:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:14
	% EndTime: 2019-10-10 00:15:15
	% DurationCPUTime: 1.03s
	% Computational Cost: add. (2921->86), mult. (2191->195), div. (456->12), fcn. (2616->9), ass. (0->91)
	t119 = qJ(1) + pkin(9);
	t115 = sin(t119);
	t109 = t115 ^ 2;
	t118 = qJ(3) + pkin(10);
	t114 = sin(t118);
	t108 = t114 ^ 2;
	t116 = cos(t118);
	t111 = 0.1e1 / t116 ^ 2;
	t161 = t108 * t111;
	t106 = t109 * t161 + 0.1e1;
	t107 = t114 * t108;
	t110 = 0.1e1 / t116;
	t160 = t110 * t114;
	t129 = qJD(3) * (t107 * t110 * t111 + t160);
	t117 = cos(t119);
	t152 = qJD(1) * t117;
	t139 = t115 * t152;
	t166 = 0.1e1 / t106 ^ 2 * (t109 * t129 + t139 * t161);
	t180 = -0.2e1 * t166;
	t104 = 0.1e1 / t106;
	t134 = 0.1e1 + t161;
	t176 = t115 * t134;
	t87 = t104 * t176;
	t179 = t115 * t87 - 0.1e1;
	t121 = cos(pkin(11));
	t151 = qJD(3) * t114;
	t136 = t121 * t151;
	t120 = sin(pkin(11));
	t155 = t117 * t120;
	t156 = t115 * t121;
	t98 = -t116 * t156 + t155;
	t92 = t98 * qJD(1) - t117 * t136;
	t154 = t117 * t121;
	t157 = t115 * t120;
	t100 = t116 * t154 + t157;
	t95 = 0.1e1 / t100 ^ 2;
	t178 = t92 * t95;
	t99 = t116 * t155 - t156;
	t168 = t95 * t99;
	t93 = t99 ^ 2;
	t90 = t93 * t95 + 0.1e1;
	t88 = 0.1e1 / t90;
	t94 = 0.1e1 / t100;
	t177 = (-t120 * t94 + t121 * t168) * t88;
	t158 = t115 * t114;
	t103 = atan2(-t158, -t116);
	t102 = cos(t103);
	t101 = sin(t103);
	t142 = t101 * t158;
	t85 = -t102 * t116 - t142;
	t82 = 0.1e1 / t85;
	t83 = 0.1e1 / t85 ^ 2;
	t175 = t104 - 0.1e1;
	t113 = t117 ^ 2;
	t149 = qJD(3) * t116;
	t143 = t83 * t149;
	t140 = t114 * t152;
	t150 = qJD(3) * t115;
	t162 = t102 * t114;
	t138 = t111 * t150;
	t78 = (-(-t115 * t149 - t140) * t110 + t108 * t138) * t104;
	t73 = (-t115 * t78 + qJD(3)) * t162 + (-t140 + (t78 - t150) * t116) * t101;
	t173 = t73 * t82 * t83;
	t81 = t108 * t113 * t83 + 0.1e1;
	t174 = (t113 * t114 * t143 + (-t113 * t173 - t83 * t139) * t108) / t81 ^ 2;
	t169 = t94 * t178;
	t137 = t120 * t151;
	t97 = -t116 * t157 - t154;
	t91 = t97 * qJD(1) - t117 * t137;
	t172 = (t91 * t168 - t93 * t169) / t90 ^ 2;
	t79 = 0.1e1 / t81;
	t171 = t79 * t83;
	t170 = t82 * t79;
	t164 = t117 * t83;
	t163 = qJD(3) * t87;
	t159 = t114 * t117;
	t153 = qJD(1) * t115;
	t148 = qJD(3) * t117;
	t147 = 0.2e1 * t173;
	t146 = 0.2e1 * t172;
	t145 = t82 * t174;
	t144 = t99 * t169;
	t141 = t104 * t108 * t110;
	t135 = 0.2e1 * t83 * t174;
	t133 = t110 * t180;
	t132 = t115 * t141;
	t131 = t134 * t117;
	t77 = (t175 * t114 * t101 - t102 * t132) * t117;
	t75 = -t179 * t162 + (-t115 + t87) * t116 * t101;
	t74 = t176 * t180 + (qJD(1) * t131 + 0.2e1 * t115 * t129) * t104;
	t1 = [t133 * t159 + (qJD(3) * t131 - t153 * t160) * t104, 0, t74, 0, 0, 0; (-t149 * t170 + (0.2e1 * t145 + (qJD(1) * t77 + t73) * t171) * t114) * t115 + (t77 * t135 * t114 + (-t77 * t143 + (t77 * t147 + ((0.2e1 * t114 * t166 - t78 * t132 - t175 * t149) * t101 + (t108 * t115 * t133 + t114 * t78 + (t107 * t138 - (t78 - 0.2e1 * t150) * t114) * t104) * t102) * t164) * t114 + (-t82 + (-(t109 - t113) * t102 * t141 + t175 * t142) * t83) * t114 * qJD(1)) * t79) * t117, 0, (-t153 * t170 + (-0.2e1 * t145 + (-qJD(3) * t75 - t73) * t171) * t117) * t116 + (t75 * t117 * t135 + (-t82 * t148 - ((-t115 * t74 - t152 * t87) * t102 + (t179 * t78 + t150 - t163) * t101) * t83 * t159 + (t117 * t147 + t83 * t153) * t75 - ((t74 - t152) * t101 + (t78 * t87 + qJD(3) + (-t78 - t163) * t115) * t102) * t116 * t164) * t79) * t114, 0, 0, 0; (t98 * t168 - t94 * t97) * t146 + ((-t99 * qJD(1) + t115 * t137) * t94 + 0.2e1 * t98 * t144 + (-t97 * t92 - (-t100 * qJD(1) + t115 * t136) * t99 - t98 * t91) * t95) * t88, 0, t116 * t148 * t177 + (-t153 * t177 + ((t94 * t146 + t178 * t88) * t120 + (-0.2e1 * t168 * t172 + (t91 * t95 - 0.2e1 * t144) * t88) * t121) * t117) * t114, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:14
	% DurationCPUTime: 1.13s
	% Computational Cost: add. (3592->98), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->94)
	t153 = qJ(3) + pkin(10);
	t147 = sin(t153);
	t140 = t147 ^ 2;
	t150 = cos(t153);
	t143 = 0.1e1 / t150 ^ 2;
	t199 = t140 * t143;
	t154 = qJ(1) + pkin(9);
	t148 = sin(t154);
	t216 = 0.2e1 * t148;
	t215 = t147 * t199;
	t152 = pkin(11) + qJ(6);
	t146 = sin(t152);
	t149 = cos(t152);
	t151 = cos(t154);
	t190 = t150 * t151;
	t129 = t146 * t148 + t149 * t190;
	t193 = t148 * t147;
	t132 = atan2(-t193, -t150);
	t131 = cos(t132);
	t130 = sin(t132);
	t179 = t130 * t193;
	t116 = -t131 * t150 - t179;
	t113 = 0.1e1 / t116;
	t123 = 0.1e1 / t129;
	t142 = 0.1e1 / t150;
	t114 = 0.1e1 / t116 ^ 2;
	t124 = 0.1e1 / t129 ^ 2;
	t214 = -0.2e1 * t147;
	t141 = t148 ^ 2;
	t135 = t141 * t199 + 0.1e1;
	t133 = 0.1e1 / t135;
	t213 = t133 - 0.1e1;
	t188 = qJD(1) * t151;
	t176 = t147 * t188;
	t186 = qJD(3) * t150;
	t187 = qJD(3) * t148;
	t107 = (-(-t148 * t186 - t176) * t142 + t187 * t199) * t133;
	t201 = t131 * t147;
	t102 = (-t107 * t148 + qJD(3)) * t201 + (-t176 + (t107 - t187) * t150) * t130;
	t212 = t102 * t113 * t114;
	t191 = t148 * t150;
	t163 = t146 * t191 + t149 * t151;
	t185 = qJD(3) * t151;
	t175 = t147 * t185;
	t108 = t163 * qJD(1) - qJD(6) * t129 + t146 * t175;
	t192 = t148 * t149;
	t128 = t146 * t190 - t192;
	t122 = t128 ^ 2;
	t121 = t122 * t124 + 0.1e1;
	t203 = t124 * t128;
	t169 = -qJD(1) * t150 + qJD(6);
	t170 = qJD(6) * t150 - qJD(1);
	t195 = t146 * t151;
	t109 = -t170 * t195 + (t169 * t148 - t175) * t149;
	t208 = t109 * t123 * t124;
	t211 = (-t108 * t203 - t122 * t208) / t121 ^ 2;
	t210 = t107 * t130;
	t209 = t107 * t147;
	t207 = t114 * t147;
	t197 = t142 * t147;
	t162 = qJD(3) * (t142 * t215 + t197);
	t167 = t140 * t148 * t188;
	t206 = (t141 * t162 + t143 * t167) / t135 ^ 2;
	t174 = 0.1e1 + t199;
	t120 = t174 * t148 * t133;
	t205 = t120 * t148;
	t204 = t123 * t146;
	t202 = t128 * t149;
	t200 = t140 * t142;
	t145 = t151 ^ 2;
	t198 = t140 * t145;
	t194 = t147 * t151;
	t189 = qJD(1) * t148;
	t112 = t114 * t198 + 0.1e1;
	t184 = 0.2e1 * (-t198 * t212 + (t145 * t147 * t186 - t167) * t114) / t112 ^ 2;
	t183 = 0.2e1 * t212;
	t182 = -0.2e1 * t211;
	t181 = t128 * t208;
	t180 = t114 * t194;
	t178 = t133 * t200;
	t173 = t147 * t184;
	t172 = t206 * t214;
	t171 = t206 * t216;
	t168 = t148 * t178;
	t166 = t174 * t151;
	t165 = t169 * t151;
	t164 = t124 * t202 - t204;
	t127 = -t149 * t191 + t195;
	t118 = 0.1e1 / t121;
	t110 = 0.1e1 / t112;
	t106 = (t213 * t147 * t130 - t131 * t168) * t151;
	t105 = -t130 * t191 + t201 + (t130 * t150 - t131 * t193) * t120;
	t103 = -t174 * t171 + (qJD(1) * t166 + t162 * t216) * t133;
	t1 = [t142 * t151 * t172 + (qJD(3) * t166 - t189 * t197) * t133, 0, t103, 0, 0, 0; (t113 * t173 + (-t113 * t186 + (qJD(1) * t106 + t102) * t207) * t110) * t148 + (t114 * t173 * t106 + (-((t107 * t168 + t213 * t186 + t172) * t130 + (t171 * t200 - t209 + (t209 + (t214 - t215) * t187) * t133) * t131) * t180 + (-t114 * t186 + t147 * t183) * t106 + (-t113 + ((-t141 + t145) * t131 * t178 + t213 * t179) * t114) * t147 * qJD(1)) * t110) * t151, 0, (t105 * t207 - t113 * t150) * t151 * t184 + ((-t113 * t189 + (-qJD(3) * t105 - t102) * t151 * t114) * t150 + (-t113 * t185 - (-t103 * t131 * t148 + t130 * t187 + t205 * t210 - t210 + (-qJD(3) * t130 - t131 * t188) * t120) * t180 + (t114 * t189 + t151 * t183) * t105 - ((t103 - t188) * t130 + ((0.1e1 - t205) * qJD(3) + (t120 - t148) * t107) * t131) * t114 * t190) * t147) * t110, 0, 0, 0; 0.2e1 * (t123 * t163 + t127 * t203) * t211 + (0.2e1 * t127 * t181 - t170 * t123 * t192 + (t147 * t187 + t165) * t204 + (t127 * t108 + t163 * t109 - t165 * t202 - (qJD(3) * t147 * t149 + t170 * t146) * t128 * t148) * t124) * t118, 0, t164 * t182 * t194 + (t164 * t150 * t185 + (-t164 * t189 + ((-qJD(6) * t123 - 0.2e1 * t181) * t149 + (-t108 * t149 + (-qJD(6) * t128 + t109) * t146) * t124) * t151) * t147) * t118, 0, 0, t182 + 0.2e1 * (-t108 * t118 * t124 + (-t118 * t208 - t124 * t211) * t128) * t128;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end