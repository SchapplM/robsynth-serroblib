% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR6
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
%   Wie in S6RPPRPR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:42
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:45
	% EndTime: 2019-10-09 23:42:45
	% DurationCPUTime: 0.76s
	% Computational Cost: add. (396->67), mult. (1839->157), div. (436->14), fcn. (2165->7), ass. (0->74)
	t88 = cos(qJ(1));
	t115 = qJD(1) * t88;
	t86 = sin(qJ(1));
	t79 = 0.1e1 / t86 ^ 2;
	t84 = t88 ^ 2;
	t122 = t79 * t84;
	t85 = sin(qJ(4));
	t73 = t85 ^ 2;
	t72 = t73 * t122 + 0.1e1;
	t77 = t86 ^ 2;
	t78 = 0.1e1 / t86;
	t95 = (-0.1e1 - 0.1e1 / t77 * t84) * t78 * t115;
	t114 = qJD(4) * t85;
	t87 = cos(qJ(4));
	t99 = t84 * t87 * t114;
	t137 = (t73 * t95 + t79 * t99) / t72 ^ 2;
	t75 = 0.1e1 / t85 ^ 2;
	t82 = t87 ^ 2;
	t123 = t75 * t82;
	t105 = 0.1e1 + t123;
	t135 = t105 * t86;
	t113 = qJD(4) * t86;
	t106 = t87 * t115;
	t71 = t77 * t123 + 0.1e1;
	t66 = 0.1e1 / t71;
	t74 = 0.1e1 / t85;
	t53 = ((-t85 * t113 + t106) * t74 - t113 * t123) * t66;
	t134 = -t53 - t113;
	t119 = t86 * t87;
	t70 = atan2(t119, t85);
	t64 = sin(t70);
	t108 = t64 * t119;
	t65 = cos(t70);
	t60 = t65 * t85 + t108;
	t57 = 0.1e1 / t60;
	t58 = 0.1e1 / t60 ^ 2;
	t133 = -0.2e1 * t87;
	t132 = t66 - 0.1e1;
	t120 = t82 * t86;
	t100 = t115 * t120;
	t121 = t82 * t84;
	t124 = t65 * t87;
	t49 = (t53 * t86 + qJD(4)) * t124 + (t134 * t85 + t106) * t64;
	t130 = t49 * t57 * t58;
	t56 = t58 * t121 + 0.1e1;
	t131 = (-t121 * t130 + (-t99 - t100) * t58) / t56 ^ 2;
	t129 = t53 * t87;
	t128 = t58 * t87;
	t127 = t58 * t88;
	t81 = t87 * t82;
	t96 = (t87 + 0.1e1 / t73 * t81) * t74;
	t126 = (-qJD(4) * t77 * t96 + t100 * t75) / t71 ^ 2;
	t125 = t64 * t86;
	t118 = t87 * t88;
	t117 = qJD(1) * t86;
	t116 = qJD(1) * t87;
	t112 = qJD(4) * t88;
	t111 = -0.2e1 * t130;
	t110 = 0.2e1 * t126;
	t109 = t58 * t118;
	t107 = t66 * t74 * t82;
	t104 = 0.1e1 + t122;
	t103 = t131 * t133;
	t102 = -0.2e1 * t74 * t126;
	t101 = t86 * t107;
	t98 = t105 * t88;
	t97 = t104 * t87;
	t68 = 0.1e1 / t72;
	t63 = t66 * t135;
	t54 = 0.1e1 / t56;
	t52 = (-t132 * t87 * t64 + t65 * t101) * t88;
	t51 = -t85 * t125 + t124 - (t65 * t119 - t64 * t85) * t63;
	t50 = t110 * t135 + (-qJD(1) * t98 + 0.2e1 * t113 * t96) * t66;
	t1 = [t102 * t118 + (-t74 * t86 * t116 - qJD(4) * t98) * t66, 0, 0, t50, 0, 0; (t57 * t103 + (-t57 * t114 + (-qJD(1) * t52 - t49) * t128) * t54) * t86 + (t58 * t103 * t52 + (((-t53 * t101 + t87 * t110 + t132 * t114) * t64 + (t102 * t120 + t129 + (-t129 + (-t75 * t81 + t133) * t113) * t66) * t65) * t109 + (t87 * t111 - t58 * t114) * t52 + (t57 + ((-t77 + t84) * t65 * t107 + t132 * t108) * t58) * t116) * t54) * t88, 0, 0, 0.2e1 * (-t51 * t128 - t57 * t85) * t88 * t131 + ((-t57 * t117 + (-qJD(4) * t51 - t49) * t127) * t85 + (t57 * t112 + (t50 * t65 * t86 + t134 * t64 - (-qJD(4) * t64 + t115 * t65 - t125 * t53) * t63) * t109 + (t88 * t111 - t58 * t117) * t51 + ((-t50 - t115) * t64 + ((t63 * t86 - 0.1e1) * qJD(4) + (t63 - t86) * t53) * t65) * t85 * t127) * t87) * t54, 0, 0; -0.2e1 * t104 * t85 * t137 + (qJD(4) * t97 + 0.2e1 * t85 * t95) * t68, 0, 0, 0.2e1 * t78 * t118 * t137 + (t78 * t85 * t112 + qJD(1) * t97) * t68, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:45
	% EndTime: 2019-10-09 23:42:46
	% DurationCPUTime: 1.11s
	% Computational Cost: add. (813->89), mult. (2519->205), div. (480->12), fcn. (2968->9), ass. (0->90)
	t126 = sin(qJ(6));
	t130 = cos(qJ(4));
	t131 = cos(qJ(1));
	t169 = t130 * t131;
	t128 = sin(qJ(1));
	t129 = cos(qJ(6));
	t172 = t128 * t129;
	t109 = -t126 * t169 - t172;
	t104 = 0.1e1 / t109;
	t105 = 0.1e1 / t109 ^ 2;
	t107 = t126 * t128 - t129 * t169;
	t178 = t107 * t126;
	t140 = -t104 * t129 + t105 * t178;
	t103 = t107 ^ 2;
	t102 = t103 * t105 + 0.1e1;
	t99 = 0.1e1 / t102;
	t191 = t140 * t99;
	t166 = qJD(4) * t128;
	t121 = t128 ^ 2;
	t127 = sin(qJ(4));
	t120 = t127 ^ 2;
	t123 = 0.1e1 / t130 ^ 2;
	t175 = t120 * t123;
	t117 = t121 * t175 + 0.1e1;
	t115 = 0.1e1 / t117;
	t122 = 0.1e1 / t130;
	t154 = t123 * t166;
	t167 = qJD(1) * t131;
	t156 = t127 * t167;
	t165 = qJD(4) * t130;
	t88 = ((-t128 * t165 - t156) * t122 - t120 * t154) * t115;
	t147 = -t88 - t166;
	t190 = t128 * t88 + qJD(4);
	t173 = t128 * t127;
	t114 = atan2(-t173, t130);
	t113 = cos(t114);
	t112 = sin(t114);
	t158 = t112 * t173;
	t97 = t113 * t130 - t158;
	t94 = 0.1e1 / t97;
	t95 = 0.1e1 / t97 ^ 2;
	t189 = t115 - 0.1e1;
	t125 = t131 ^ 2;
	t155 = t128 * t167;
	t159 = t95 * t165;
	t177 = t113 * t127;
	t83 = -t190 * t177 + (t130 * t147 - t156) * t112;
	t187 = t83 * t94 * t95;
	t93 = t120 * t125 * t95 + 0.1e1;
	t188 = (t125 * t127 * t159 + (-t125 * t187 - t155 * t95) * t120) / t93 ^ 2;
	t91 = 0.1e1 / t93;
	t186 = t91 * t95;
	t185 = t94 * t91;
	t179 = t105 * t107;
	t145 = qJD(1) * t130 + qJD(6);
	t164 = qJD(4) * t131;
	t139 = t127 * t164 + t128 * t145;
	t146 = qJD(6) * t130 + qJD(1);
	t170 = t129 * t131;
	t90 = t126 * t139 - t146 * t170;
	t183 = t104 * t105 * t90;
	t141 = t146 * t126;
	t89 = t129 * t139 + t131 * t141;
	t184 = 0.1e1 / t102 ^ 2 * (-t103 * t183 + t179 * t89);
	t181 = t131 * t95;
	t176 = t120 * t122;
	t174 = t122 * t127;
	t171 = t128 * t130;
	t168 = qJD(1) * t128;
	t163 = 0.2e1 * t187;
	t162 = -0.2e1 * t184;
	t119 = t127 * t120;
	t142 = t119 * t122 * t123 + t174;
	t161 = 0.2e1 / t117 ^ 2 * (qJD(4) * t121 * t142 + t155 * t175);
	t160 = t94 * t188;
	t157 = t115 * t176;
	t153 = 0.2e1 * t95 * t188;
	t152 = 0.1e1 + t175;
	t151 = -0.2e1 * t107 * t183;
	t150 = t128 * t161;
	t149 = t127 * t161;
	t144 = t128 * t157;
	t143 = t152 * t131;
	t111 = t126 * t171 - t170;
	t110 = -t126 * t131 - t129 * t171;
	t101 = t152 * t128 * t115;
	t87 = (t112 * t127 * t189 + t113 * t144) * t131;
	t86 = -t112 * t171 - t177 - (-t112 * t130 - t113 * t173) * t101;
	t84 = t152 * t150 + (-qJD(1) * t143 - 0.2e1 * t142 * t166) * t115;
	t1 = [t122 * t131 * t149 + (-qJD(4) * t143 + t168 * t174) * t115, 0, 0, t84, 0, 0; (-t165 * t185 + (0.2e1 * t160 + (qJD(1) * t87 + t83) * t186) * t127) * t128 + (t87 * t153 * t127 + (-t87 * t159 + (t87 * t163 + ((t144 * t88 - t165 * t189 + t149) * t112 + (t150 * t176 + t88 * t127 + (-t119 * t154 - (t88 + 0.2e1 * t166) * t127) * t115) * t113) * t181) * t127 + (-t94 + (-(-t121 + t125) * t113 * t157 + t189 * t158) * t95) * t127 * qJD(1)) * t91) * t131, 0, 0, (-t168 * t185 + (-0.2e1 * t160 + (-qJD(4) * t86 - t83) * t186) * t131) * t130 + (t86 * t131 * t153 + (-t94 * t164 - ((t101 * t167 - t128 * t84) * t113 + (-t190 * t101 - t147) * t112) * t127 * t181 + (t131 * t163 + t168 * t95) * t86) * t91 - ((-t84 - t167) * t112 + (-t101 * t147 - t190) * t113) * t169 * t186) * t127, 0, 0; 0.2e1 * (-t104 * t110 - t111 * t179) * t184 + (t111 * t151 + (-t110 * t90 + t111 * t89 + t146 * t107 * t172 + (-t127 * t166 + t131 * t145) * t178) * t105 + (-t145 * t170 + (qJD(4) * t127 * t129 + t141) * t128) * t104) * t99, 0, 0, t130 * t164 * t191 + (-t168 * t191 + (t140 * t162 + ((qJD(6) * t104 + t151) * t126 + (t126 * t89 + (qJD(6) * t107 + t90) * t129) * t105) * t99) * t131) * t127, 0, t162 + 0.2e1 * (t105 * t89 * t99 + (-t105 * t184 - t183 * t99) * t107) * t107;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end