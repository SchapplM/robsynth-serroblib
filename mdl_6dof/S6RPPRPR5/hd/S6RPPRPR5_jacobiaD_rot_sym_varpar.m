% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR5
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
%   Wie in S6RPPRPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:05
	% EndTime: 2019-10-09 23:41:05
	% DurationCPUTime: 0.92s
	% Computational Cost: add. (514->81), mult. (2191->187), div. (456->12), fcn. (2616->9), ass. (0->88)
	t102 = sin(qJ(4));
	t94 = 0.1e1 / t102 ^ 2;
	t104 = cos(qJ(4));
	t98 = t104 ^ 2;
	t151 = t94 * t98;
	t103 = sin(qJ(1));
	t126 = 0.1e1 + t151;
	t96 = t103 ^ 2;
	t92 = t96 * t151 + 0.1e1;
	t89 = 0.1e1 / t92;
	t115 = t126 * t89;
	t76 = t103 * t115;
	t163 = t103 * t76 - 0.1e1;
	t101 = cos(pkin(9));
	t105 = cos(qJ(1));
	t133 = qJD(4) * t105;
	t122 = t104 * t133;
	t100 = sin(pkin(9));
	t141 = t105 * t100;
	t143 = t103 * t101;
	t84 = -t102 * t143 - t141;
	t78 = t84 * qJD(1) + t101 * t122;
	t140 = t105 * t101;
	t144 = t103 * t100;
	t86 = t102 * t140 - t144;
	t81 = 0.1e1 / t86 ^ 2;
	t162 = t78 * t81;
	t161 = t104 * t151;
	t142 = t103 * t104;
	t91 = atan2(t142, t102);
	t87 = sin(t91);
	t127 = t87 * t142;
	t88 = cos(t91);
	t74 = t88 * t102 + t127;
	t70 = 0.1e1 / t74;
	t80 = 0.1e1 / t86;
	t93 = 0.1e1 / t102;
	t71 = 0.1e1 / t74 ^ 2;
	t160 = t89 - 0.1e1;
	t135 = qJD(4) * t103;
	t137 = qJD(1) * t105;
	t149 = t102 * t87;
	t64 = ((-t102 * t135 + t104 * t137) * t93 - t135 * t151) * t89;
	t59 = (-t64 - t135) * t149 + (t87 * t137 + (t103 * t64 + qJD(4)) * t88) * t104;
	t159 = t59 * t70 * t71;
	t85 = t102 * t141 + t143;
	t153 = t81 * t85;
	t154 = t80 * t162;
	t79 = t85 ^ 2;
	t73 = t79 * t81 + 0.1e1;
	t83 = -t102 * t144 + t140;
	t77 = t83 * qJD(1) + t100 * t122;
	t158 = (t77 * t153 - t79 * t154) / t73 ^ 2;
	t99 = t105 ^ 2;
	t150 = t98 * t99;
	t67 = t71 * t150 + 0.1e1;
	t65 = 0.1e1 / t67;
	t156 = t65 * t71;
	t114 = (t104 + t161) * t93;
	t116 = t103 * t98 * t137;
	t155 = (-t114 * t96 * qJD(4) + t94 * t116) / t92 ^ 2;
	t152 = t89 * t93;
	t147 = t105 * t71;
	t146 = qJD(4) * t76;
	t145 = qJD(4) * t89;
	t139 = qJD(1) * t103;
	t138 = qJD(1) * t104;
	t136 = qJD(4) * t102;
	t134 = qJD(4) * t104;
	t132 = -0.2e1 * (-t150 * t159 + (-t102 * t99 * t134 - t116) * t71) / t67 ^ 2;
	t131 = -0.2e1 * t159;
	t130 = 0.2e1 * t158;
	t129 = 0.2e1 * t155;
	t128 = t103 * t152;
	t125 = t160 * t104;
	t124 = t65 * t136;
	t123 = t103 * t134;
	t121 = t70 * t132;
	t120 = t71 * t132;
	t119 = -0.2e1 * t93 * t155;
	t118 = 0.2e1 * t85 * t154;
	t117 = t98 * t128;
	t68 = 0.1e1 / t73;
	t113 = (-t100 * t80 + t101 * t153) * t68;
	t63 = (t88 * t117 - t87 * t125) * t105;
	t61 = -t163 * t88 * t104 + (-t103 + t76) * t149;
	t60 = -t115 * t137 + (0.2e1 * t114 * t145 + t126 * t129) * t103;
	t1 = [-t128 * t138 + (-qJD(4) * t115 + t104 * t119) * t105, 0, 0, t60, 0, 0; (-t70 * t124 + (t121 + (-qJD(1) * t63 - t59) * t156) * t104) * t103 + (-t63 * t71 * t124 + (t63 * t120 + (t63 * t131 + ((t104 * t129 - t64 * t117 + t160 * t136) * t87 + (-t64 * t125 + (t98 * t119 + (-0.2e1 * t104 - t161) * t145) * t103) * t88) * t147) * t65) * t104 + (t70 + ((-t96 + t99) * t98 * t88 * t152 + t160 * t127) * t71) * t65 * t138) * t105, 0, 0, (-t70 * t65 * t139 + (t121 + (-qJD(4) * t61 - t59) * t156) * t105) * t102 + (t61 * t105 * t120 + (t70 * t133 + (t105 * t131 - t71 * t139) * t61 + (((t103 * t60 - t137 * t76) * t88 + (t163 * t64 - t135 + t146) * t87) * t104 + ((-t60 - t137) * t87 + (t64 * t76 - qJD(4) + (-t64 + t146) * t103) * t88) * t102) * t147) * t65) * t104, 0, 0; (t84 * t153 - t80 * t83) * t130 + ((-t85 * qJD(1) - t100 * t123) * t80 + t84 * t118 + (-t83 * t78 - (-t86 * qJD(1) - t101 * t123) * t85 - t84 * t77) * t81) * t68, 0, 0, t102 * t113 * t133 + (t113 * t139 + ((-0.2e1 * t80 * t158 - t68 * t162) * t100 + (t130 * t153 + (-t77 * t81 + t118) * t68) * t101) * t105) * t104, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:05
	% EndTime: 2019-10-09 23:41:05
	% DurationCPUTime: 0.95s
	% Computational Cost: add. (972->92), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->95)
	t135 = sin(qJ(1));
	t134 = sin(qJ(4));
	t128 = 0.1e1 / t134 ^ 2;
	t136 = cos(qJ(4));
	t132 = t136 ^ 2;
	t185 = t128 * t132;
	t159 = 0.1e1 + t185;
	t198 = t135 * t159;
	t179 = t135 * t136;
	t121 = atan2(t179, t134);
	t118 = cos(t121);
	t117 = sin(t121);
	t165 = t117 * t179;
	t104 = t118 * t134 + t165;
	t101 = 0.1e1 / t104;
	t126 = pkin(9) + qJ(6);
	t125 = cos(t126);
	t137 = cos(qJ(1));
	t182 = t134 * t137;
	t163 = t125 * t182;
	t124 = sin(t126);
	t181 = t135 * t124;
	t114 = t163 - t181;
	t108 = 0.1e1 / t114;
	t127 = 0.1e1 / t134;
	t102 = 0.1e1 / t104 ^ 2;
	t109 = 0.1e1 / t114 ^ 2;
	t130 = t135 ^ 2;
	t122 = t130 * t185 + 0.1e1;
	t119 = 0.1e1 / t122;
	t197 = t119 - 0.1e1;
	t180 = t135 * t125;
	t113 = t124 * t182 + t180;
	t107 = t113 ^ 2;
	t190 = t109 * t113;
	t152 = qJD(1) * t134 + qJD(6);
	t153 = qJD(6) * t134 + qJD(1);
	t171 = qJD(4) * t137;
	t160 = t136 * t171;
	t186 = t124 * t137;
	t93 = -t153 * t186 + (-t152 * t135 + t160) * t125;
	t193 = t108 * t109 * t93;
	t175 = qJD(1) * t137;
	t92 = -qJD(6) * t163 - t124 * t160 - t125 * t175 + t152 * t181;
	t97 = t107 * t109 + 0.1e1;
	t196 = (-t107 * t193 - t92 * t190) / t97 ^ 2;
	t133 = t137 ^ 2;
	t184 = t132 * t133;
	t100 = t102 * t184 + 0.1e1;
	t98 = 0.1e1 / t100;
	t195 = t102 * t98;
	t173 = qJD(4) * t135;
	t161 = t128 * t173;
	t162 = t136 * t175;
	t94 = ((-t134 * t173 + t162) * t127 - t132 * t161) * t119;
	t154 = -t94 - t173;
	t155 = t135 * t94 + qJD(4);
	t187 = t118 * t136;
	t88 = t155 * t187 + (t154 * t134 + t162) * t117;
	t194 = t101 * t102 * t88;
	t131 = t136 * t132;
	t148 = (t128 * t131 + t136) * t127;
	t183 = t132 * t135;
	t150 = t175 * t183;
	t192 = (-t148 * t130 * qJD(4) + t128 * t150) / t122 ^ 2;
	t191 = t108 * t124;
	t189 = t113 * t125;
	t188 = t117 * t134;
	t178 = t136 * t137;
	t177 = qJD(1) * t135;
	t176 = qJD(1) * t136;
	t174 = qJD(4) * t134;
	t172 = qJD(4) * t136;
	t170 = -0.2e1 * (-t184 * t194 + (-t133 * t134 * t172 - t150) * t102) / t100 ^ 2;
	t169 = 0.2e1 * t196;
	t168 = -0.2e1 * t194;
	t167 = 0.2e1 * t192;
	t166 = t98 * t174;
	t164 = t119 * t127 * t132;
	t158 = t101 * t170;
	t157 = 0.2e1 * t113 * t193;
	t156 = -0.2e1 * t127 * t192;
	t151 = t135 * t164;
	t149 = t159 * t137;
	t147 = t109 * t189 - t191;
	t95 = 0.1e1 / t97;
	t146 = t147 * t95;
	t145 = -t135 * t172 - t152 * t137;
	t112 = -t134 * t180 - t186;
	t111 = t125 * t137 - t134 * t181;
	t106 = t119 * t198;
	t91 = (-t197 * t136 * t117 + t118 * t151) * t137;
	t90 = -t135 * t188 + t187 - (t118 * t179 - t188) * t106;
	t89 = t167 * t198 + (-qJD(1) * t149 + 0.2e1 * t148 * t173) * t119;
	t1 = [t156 * t178 + (-t127 * t135 * t176 - qJD(4) * t149) * t119, 0, 0, t89, 0, 0; (-t101 * t166 + (t158 + (-qJD(1) * t91 - t88) * t195) * t136) * t135 + (t91 * t136 * t98 * t168 + (-t91 * t166 + (t91 * t170 + ((t136 * t167 - t94 * t151 + t197 * t174) * t117 + (t156 * t183 + t136 * t94 + (-t131 * t161 + (-t94 - 0.2e1 * t173) * t136) * t119) * t118) * t98 * t137) * t136) * t102 + (t101 + ((-t130 + t133) * t118 * t164 + t197 * t165) * t102) * t98 * t176) * t137, 0, 0, (-t101 * t98 * t177 + (t158 + (-qJD(4) * t90 - t88) * t195) * t137) * t134 + (t90 * t137 * t102 * t170 + ((qJD(4) * t101 + t90 * t168) * t137 + (-t90 * t177 + ((-t106 * t175 + t135 * t89) * t118 + ((t106 * t135 - 0.1e1) * t94 + (t106 - t135) * qJD(4)) * t117) * t178) * t102) * t98 + ((-t89 - t175) * t117 + (-t154 * t106 - t155) * t118) * t182 * t195) * t136, 0, 0; (-t108 * t111 + t112 * t190) * t169 + (t112 * t157 - t153 * t108 * t180 + t145 * t191 + (-t153 * t113 * t181 - t111 * t93 + t112 * t92 - t145 * t189) * t109) * t95, 0, 0, t134 * t146 * t171 + (t146 * t177 + (t147 * t169 + ((qJD(6) * t108 + t157) * t125 + (t125 * t92 + (qJD(6) * t113 - t93) * t124) * t109) * t95) * t137) * t136, 0, -0.2e1 * t196 + 0.2e1 * (-t109 * t92 * t95 + (-t109 * t196 - t95 * t193) * t113) * t113;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end