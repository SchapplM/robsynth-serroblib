% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR8
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
%   Wie in S6RPRPPR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:27
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR8_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR8_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:24
	% DurationCPUTime: 0.76s
	% Computational Cost: add. (587->69), mult. (1835->162), div. (470->13), fcn. (2177->7), ass. (0->72)
	t84 = cos(qJ(1));
	t111 = qJD(1) * t84;
	t81 = sin(qJ(3));
	t70 = 0.1e1 / t81;
	t133 = 0.2e1 * t70;
	t132 = 0.2e1 * t84;
	t108 = qJD(3) * t84;
	t83 = cos(qJ(3));
	t112 = qJD(1) * t83;
	t82 = sin(qJ(1));
	t102 = t82 * t112;
	t71 = 0.1e1 / t81 ^ 2;
	t78 = t83 ^ 2;
	t117 = t71 * t78;
	t80 = t84 ^ 2;
	t69 = t80 * t117 + 0.1e1;
	t66 = 0.1e1 / t69;
	t50 = ((t81 * t108 + t102) * t70 + t108 * t117) * t66;
	t130 = -t50 + t108;
	t114 = t84 * t83;
	t63 = atan2(-t114, t81);
	t61 = sin(t63);
	t62 = cos(t63);
	t57 = -t61 * t114 + t62 * t81;
	t54 = 0.1e1 / t57;
	t74 = 0.1e1 / t82;
	t55 = 0.1e1 / t57 ^ 2;
	t129 = t66 - 0.1e1;
	t109 = qJD(3) * t83;
	t73 = t82 ^ 2;
	t116 = t73 * t78;
	t120 = t62 * t83;
	t46 = (-t50 * t84 + qJD(3)) * t120 + (t130 * t81 + t102) * t61;
	t127 = t46 * t54 * t55;
	t53 = t55 * t116 + 0.1e1;
	t96 = t78 * t82 * t111;
	t128 = (-t116 * t127 + (-t73 * t81 * t109 + t96) * t55) / t53 ^ 2;
	t126 = t50 * t83;
	t125 = t55 * t82;
	t124 = t55 * t83;
	t118 = t70 * t83;
	t72 = t70 * t71;
	t77 = t83 * t78;
	t92 = qJD(3) * (-t72 * t77 - t118);
	t123 = (-t71 * t96 + t80 * t92) / t69 ^ 2;
	t115 = 0.1e1 / t82 ^ 2 * t80;
	t68 = t71 * t115 + 0.1e1;
	t93 = (-0.1e1 - 0.1e1 / t73 * t80) * t74 * t111;
	t122 = (-t72 * t109 * t115 + t71 * t93) / t68 ^ 2;
	t121 = t61 * t84;
	t119 = t70 * t78;
	t113 = qJD(1) * t82;
	t110 = qJD(3) * t81;
	t107 = -0.2e1 * t127;
	t106 = t83 * t128;
	t105 = t82 * t124;
	t104 = t83 * t123;
	t103 = t66 * t119;
	t101 = 0.1e1 + t117;
	t100 = -0.1e1 - t115;
	t99 = t123 * t132;
	t98 = t84 * t103;
	t97 = t129 * t83 * t61;
	t95 = t101 * t82;
	t94 = t100 * t83 * t71;
	t64 = 0.1e1 / t68;
	t60 = t101 * t84 * t66;
	t51 = 0.1e1 / t53;
	t49 = (-t62 * t98 - t97) * t82;
	t48 = t81 * t121 + t120 + (-t62 * t114 - t61 * t81) * t60;
	t47 = -t101 * t99 + (-qJD(1) * t95 + t92 * t132) * t66;
	t1 = [-0.2e1 * t82 * t70 * t104 + (-qJD(3) * t95 + t111 * t118) * t66, 0, t47, 0, 0, 0; (0.2e1 * t54 * t106 + (t54 * t110 + (qJD(1) * t49 + t46) * t124) * t51) * t84 + (-0.2e1 * t55 * t106 * t49 + (((t129 * t110 + t50 * t98 + 0.2e1 * t104) * t61 + (t99 * t119 + t126 + (-t126 + (t71 * t77 + 0.2e1 * t83) * t108) * t66) * t62) * t105 + (t83 * t107 - t55 * t110) * t49 + (t54 + ((t73 - t80) * t62 * t103 - t84 * t97) * t55) * t112) * t51) * t82, 0, 0.2e1 * (-t48 * t124 - t54 * t81) * t82 * t128 + ((t54 * t111 + (-qJD(3) * t48 - t46) * t125) * t81 + (t82 * qJD(3) * t54 + (-t47 * t62 * t84 + t130 * t61 + (-qJD(3) * t61 + t113 * t62 + t121 * t50) * t60) * t105 + (t82 * t107 + t55 * t111) * t48 + ((-t47 - t113) * t61 + ((t60 * t84 - 0.1e1) * qJD(3) + (-t60 + t84) * t50) * t62) * t81 * t125) * t83) * t51, 0, 0, 0; t100 * t122 * t133 + (qJD(3) * t94 + t93 * t133) * t64, 0, -0.2e1 * t74 * t71 * t114 * t122 + ((-0.2e1 * t72 * t78 - t70) * t74 * t108 + qJD(1) * t94) * t64, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:24
	% DurationCPUTime: 1.00s
	% Computational Cost: add. (624->89), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->92)
	t128 = sin(qJ(6));
	t131 = cos(qJ(6));
	t133 = cos(qJ(1));
	t171 = t133 * t131;
	t130 = sin(qJ(1));
	t132 = cos(qJ(3));
	t174 = t130 * t132;
	t108 = t128 * t174 - t171;
	t194 = 0.2e1 * t108;
	t127 = t133 ^ 2;
	t129 = sin(qJ(3));
	t122 = t129 ^ 2;
	t125 = 0.1e1 / t132 ^ 2;
	t177 = t122 * t125;
	t118 = t127 * t177 + 0.1e1;
	t121 = t129 * t122;
	t124 = 0.1e1 / t132;
	t176 = t124 * t129;
	t142 = qJD(3) * (t121 * t124 * t125 + t176);
	t169 = qJD(1) * t133;
	t160 = t130 * t169;
	t184 = 0.1e1 / t118 ^ 2 * (t127 * t142 - t160 * t177);
	t193 = -0.2e1 * t184;
	t166 = qJD(3) * t133;
	t115 = 0.1e1 / t118;
	t159 = t125 * t166;
	t170 = qJD(1) * t130;
	t161 = t129 * t170;
	t89 = ((t132 * t166 - t161) * t124 + t122 * t159) * t115;
	t152 = -t89 + t166;
	t156 = 0.1e1 + t177;
	t192 = t133 * t156;
	t172 = t133 * t129;
	t117 = atan2(t172, t132);
	t113 = sin(t117);
	t114 = cos(t117);
	t98 = t113 * t172 + t114 * t132;
	t95 = 0.1e1 / t98;
	t173 = t133 * t128;
	t144 = t131 * t174 + t173;
	t105 = 0.1e1 / t144;
	t106 = 0.1e1 / t144 ^ 2;
	t96 = 0.1e1 / t98 ^ 2;
	t191 = t115 - 0.1e1;
	t123 = t130 ^ 2;
	t167 = qJD(3) * t132;
	t163 = t96 * t167;
	t153 = t133 * t89 - qJD(3);
	t179 = t114 * t129;
	t84 = t153 * t179 + (t152 * t132 - t161) * t113;
	t189 = t84 * t95 * t96;
	t94 = t123 * t122 * t96 + 0.1e1;
	t190 = (t123 * t129 * t163 + (-t123 * t189 + t96 * t160) * t122) / t94 ^ 2;
	t92 = 0.1e1 / t94;
	t188 = t92 * t96;
	t187 = t95 * t92;
	t104 = t108 ^ 2;
	t103 = t104 * t106 + 0.1e1;
	t182 = t106 * t108;
	t151 = qJD(6) * t132 + qJD(1);
	t146 = t151 * t128;
	t150 = qJD(1) * t132 + qJD(6);
	t91 = -t150 * t171 + (qJD(3) * t129 * t131 + t146) * t130;
	t185 = t105 * t106 * t91;
	t143 = t130 * t131 + t132 * t173;
	t168 = qJD(3) * t130;
	t90 = -t129 * t128 * t168 + t143 * qJD(1) + t144 * qJD(6);
	t186 = 0.1e1 / t103 ^ 2 * (t104 * t185 + t90 * t182);
	t183 = t105 * t128;
	t181 = t108 * t131;
	t180 = t113 * t133;
	t178 = t122 * t124;
	t175 = t129 * t130;
	t165 = 0.2e1 * t189;
	t164 = -0.2e1 * t186;
	t162 = t133 * t178;
	t158 = -0.2e1 * t95 * t190;
	t157 = 0.2e1 * t96 * t190;
	t155 = t185 * t194;
	t154 = 0.2e1 * t129 * t184;
	t149 = t115 * t162;
	t148 = t191 * t129 * t113;
	t147 = t156 * t130;
	t145 = t106 * t181 - t183;
	t141 = t129 * t166 + t150 * t130;
	t112 = t130 * t128 - t132 * t171;
	t102 = t115 * t192;
	t100 = 0.1e1 / t103;
	t88 = (-t114 * t149 + t148) * t130;
	t87 = t132 * t180 - t179 + (-t113 * t132 + t114 * t172) * t102;
	t85 = t192 * t193 + (-qJD(1) * t147 + 0.2e1 * t133 * t142) * t115;
	t1 = [t130 * t124 * t154 + (-qJD(3) * t147 - t169 * t176) * t115, 0, t85, 0, 0, 0; (t167 * t187 + (t158 + (-qJD(1) * t88 - t84) * t188) * t129) * t133 + (t88 * t157 * t129 + (-t88 * t163 + (t88 * t165 + ((-t89 * t149 - t191 * t167 + t154) * t113 + (t162 * t193 + t129 * t89 + (t121 * t159 - (t89 - 0.2e1 * t166) * t129) * t115) * t114) * t96 * t130) * t129 + (-t95 + (-(t123 - t127) * t115 * t114 * t178 - t133 * t148) * t96) * t129 * qJD(1)) * t92) * t130, 0, (t169 * t187 + (t158 + (-qJD(3) * t87 - t84) * t188) * t130) * t132 + (t87 * t130 * t157 + (-t95 * t168 - (t114 * t133 * t85 - t152 * t113 + (qJD(3) * t113 - t114 * t170 - t180 * t89) * t102) * t96 * t175 + (t130 * t165 - t96 * t169) * t87) * t92 - ((-t85 - t170) * t113 + (t152 * t102 + t153) * t114) * t174 * t188) * t129, 0, 0, 0; 0.2e1 * (-t105 * t143 - t112 * t182) * t186 + (t112 * t155 + t151 * t105 * t171 - t141 * t183 + (t133 * t108 * t146 + t112 * t90 + t141 * t181 + t143 * t91) * t106) * t100, 0, t145 * t164 * t175 + (t145 * t130 * t167 + (t145 * t169 + ((-qJD(6) * t105 + t155) * t131 + (t131 * t90 + (-qJD(6) * t108 - t91) * t128) * t106) * t130) * t129) * t100, 0, 0, t164 + (t90 * t106 * t100 + (t100 * t185 - t106 * t186) * t108) * t194;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end