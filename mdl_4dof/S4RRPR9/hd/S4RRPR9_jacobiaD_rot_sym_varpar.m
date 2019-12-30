% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S4RRPR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 13:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4RRPR9_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR9_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_jacobiaD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:46:42
	% EndTime: 2019-12-29 13:46:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:46:37
	% EndTime: 2019-12-29 13:46:37
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:46:44
	% EndTime: 2019-12-29 13:46:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:46:44
	% EndTime: 2019-12-29 13:46:45
	% DurationCPUTime: 1.42s
	% Computational Cost: add. (892->82), mult. (2191->191), div. (456->12), fcn. (2616->9), ass. (0->85)
	t100 = sin(qJ(2));
	t92 = t100 ^ 2;
	t102 = cos(qJ(2));
	t95 = 0.1e1 / t102 ^ 2;
	t141 = t92 * t95;
	t101 = sin(qJ(1));
	t122 = 0.1e1 + t141;
	t93 = t101 ^ 2;
	t90 = t93 * t141 + 0.1e1;
	t88 = 0.1e1 / t90;
	t113 = t122 * t88;
	t74 = t101 * t113;
	t157 = t101 * t74 - 0.1e1;
	t103 = cos(qJ(1));
	t127 = qJD(2) * t103;
	t118 = t100 * t127;
	t134 = t101 * t102;
	t98 = sin(pkin(7));
	t99 = cos(pkin(7));
	t82 = t103 * t98 - t99 * t134;
	t76 = t82 * qJD(1) - t99 * t118;
	t133 = t102 * t103;
	t84 = t101 * t98 + t99 * t133;
	t78 = 0.1e1 / t84;
	t79 = 0.1e1 / t84 ^ 2;
	t80 = t78 * t79;
	t145 = t76 * t80;
	t81 = -t103 * t99 - t98 * t134;
	t75 = t81 * qJD(1) - t98 * t118;
	t146 = t75 * t79;
	t83 = -t101 * t99 + t98 * t133;
	t77 = t83 ^ 2;
	t72 = t77 * t79 + 0.1e1;
	t156 = (-t77 * t145 + t83 * t146) / t72 ^ 2;
	t155 = t100 * t141;
	t143 = t83 * t99;
	t112 = t79 * t143 - t78 * t98;
	t70 = 0.1e1 / t72;
	t154 = t112 * t70;
	t135 = t101 * t100;
	t87 = atan2(-t135, -t102);
	t85 = sin(t87);
	t123 = t85 * t135;
	t86 = cos(t87);
	t69 = -t102 * t86 - t123;
	t66 = 0.1e1 / t69;
	t94 = 0.1e1 / t102;
	t67 = 0.1e1 / t69 ^ 2;
	t153 = 0.2e1 * t100;
	t152 = t88 - 0.1e1;
	t130 = qJD(1) * t103;
	t114 = t101 * t92 * t130;
	t128 = qJD(2) * t102;
	t97 = t103 ^ 2;
	t140 = t92 * t97;
	t129 = qJD(2) * t101;
	t138 = t102 * t85;
	t62 = (-(-t100 * t130 - t101 * t128) * t94 + t129 * t141) * t88;
	t57 = (t62 - t129) * t138 + (-t85 * t130 + (-t101 * t62 + qJD(2)) * t86) * t100;
	t150 = t57 * t66 * t67;
	t65 = t67 * t140 + 0.1e1;
	t151 = (-t140 * t150 + (t100 * t97 * t128 - t114) * t67) / t65 ^ 2;
	t63 = 0.1e1 / t65;
	t148 = t63 * t67;
	t111 = qJD(2) * (t100 + t155) * t94;
	t147 = (t111 * t93 + t95 * t114) / t90 ^ 2;
	t144 = t82 * t83;
	t142 = t88 * t94;
	t137 = t103 * t67;
	t136 = qJD(2) * t74;
	t132 = qJD(1) * t100;
	t131 = qJD(1) * t101;
	t126 = 0.2e1 * t150;
	t125 = t66 * t151;
	t124 = t101 * t142;
	t121 = t100 * t152;
	t120 = t63 * t128;
	t119 = t100 * t129;
	t117 = 0.2e1 * t67 * t151;
	t116 = -0.2e1 * t94 * t147;
	t115 = t92 * t124;
	t61 = (-t86 * t115 + t85 * t121) * t103;
	t59 = (-t101 + t74) * t138 - t157 * t86 * t100;
	t58 = t113 * t130 + 0.2e1 * (t111 * t88 - t122 * t147) * t101;
	t1 = [-t124 * t132 + (qJD(2) * t113 + t100 * t116) * t103, t58, 0, 0; (-t66 * t120 + (0.2e1 * t125 + (qJD(1) * t61 + t57) * t148) * t100) * t101 + (-t61 * t67 * t120 + (t61 * t117 + (t61 * t126 + ((-t62 * t115 - t152 * t128 + t147 * t153) * t85 + (-t62 * t121 + (t92 * t116 + (t153 + t155) * t88 * qJD(2)) * t101) * t86) * t137) * t63) * t100 + (-t66 + (-(t93 - t97) * t92 * t86 * t142 + t152 * t123) * t67) * t63 * t132) * t103, (-t66 * t63 * t131 + (-0.2e1 * t125 + (-qJD(2) * t59 - t57) * t148) * t103) * t102 + (t59 * t103 * t117 + (-t66 * t127 - ((-t101 * t58 - t130 * t74) * t86 + (t157 * t62 + t129 - t136) * t85) * t100 * t137 + (t103 * t126 + t67 * t131) * t59) * t63 - ((t58 - t130) * t85 + (t62 * t74 + qJD(2) + (-t62 - t136) * t101) * t86) * t133 * t148) * t100, 0, 0; 0.2e1 * (t79 * t144 - t78 * t81) * t156 + ((-t83 * qJD(1) + t98 * t119) * t78 + 0.2e1 * t144 * t145 + (-t81 * t76 - (-t84 * qJD(1) + t99 * t119) * t83 - t82 * t75) * t79) * t70, t102 * t127 * t154 + (-t131 * t154 + (-0.2e1 * t112 * t156 + (t99 * t146 + (-0.2e1 * t80 * t143 + t79 * t98) * t76) * t70) * t103) * t100, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:46:44
	% EndTime: 2019-12-29 13:46:45
	% DurationCPUTime: 1.52s
	% Computational Cost: add. (1350->93), mult. (2519->209), div. (480->12), fcn. (2968->9), ass. (0->94)
	t133 = sin(qJ(1));
	t127 = t133 ^ 2;
	t132 = sin(qJ(2));
	t126 = t132 ^ 2;
	t134 = cos(qJ(2));
	t129 = 0.1e1 / t134 ^ 2;
	t182 = t126 * t129;
	t119 = t127 * t182 + 0.1e1;
	t125 = t132 * t126;
	t128 = 0.1e1 / t134;
	t179 = t128 * t132;
	t144 = qJD(2) * (t125 * t128 * t129 + t179);
	t135 = cos(qJ(1));
	t171 = qJD(1) * t135;
	t180 = t126 * t133;
	t148 = t171 * t180;
	t188 = (t127 * t144 + t129 * t148) / t119 ^ 2;
	t198 = -0.2e1 * t188;
	t124 = pkin(7) + qJ(4);
	t123 = cos(t124);
	t173 = t134 * t135;
	t122 = sin(t124);
	t177 = t133 * t122;
	t112 = t123 * t173 + t177;
	t107 = 0.1e1 / t112 ^ 2;
	t176 = t133 * t123;
	t111 = t122 * t173 - t176;
	t185 = t111 * t123;
	t106 = 0.1e1 / t112;
	t187 = t106 * t122;
	t146 = t107 * t185 - t187;
	t105 = t111 ^ 2;
	t98 = t105 * t107 + 0.1e1;
	t96 = 0.1e1 / t98;
	t197 = t146 * t96;
	t155 = 0.1e1 + t182;
	t196 = t133 * t155;
	t175 = t133 * t132;
	t116 = atan2(-t175, -t134);
	t114 = cos(t116);
	t113 = sin(t116);
	t161 = t113 * t175;
	t102 = -t114 * t134 - t161;
	t99 = 0.1e1 / t102;
	t100 = 0.1e1 / t102 ^ 2;
	t117 = 0.1e1 / t119;
	t195 = t117 - 0.1e1;
	t131 = t135 ^ 2;
	t169 = qJD(2) * t134;
	t181 = t126 * t131;
	t170 = qJD(2) * t133;
	t157 = t129 * t170;
	t158 = t132 * t171;
	t92 = (-(-t133 * t169 - t158) * t128 + t126 * t157) * t117;
	t152 = t92 - t170;
	t153 = -t133 * t92 + qJD(2);
	t184 = t114 * t132;
	t86 = t153 * t184 + (t152 * t134 - t158) * t113;
	t191 = t99 * t100 * t86;
	t95 = t100 * t181 + 0.1e1;
	t194 = (-t181 * t191 + (t131 * t132 * t169 - t148) * t100) / t95 ^ 2;
	t186 = t107 * t111;
	t150 = -qJD(1) * t134 + qJD(4);
	t151 = qJD(4) * t134 - qJD(1);
	t168 = qJD(2) * t135;
	t156 = t132 * t168;
	t183 = t122 * t135;
	t91 = -t151 * t183 + (t150 * t133 - t156) * t123;
	t190 = t106 * t107 * t91;
	t174 = t133 * t134;
	t145 = t122 * t174 + t123 * t135;
	t90 = t145 * qJD(1) - qJD(4) * t112 + t122 * t156;
	t193 = (-t105 * t190 - t90 * t186) / t98 ^ 2;
	t93 = 0.1e1 / t95;
	t192 = t100 * t93;
	t178 = t132 * t135;
	t172 = qJD(1) * t133;
	t167 = 0.2e1 * t194;
	t166 = -0.2e1 * t193;
	t165 = 0.2e1 * t191;
	t164 = t99 * t194;
	t163 = t111 * t190;
	t162 = t93 * t169;
	t160 = t117 * t126 * t128;
	t154 = t128 * t198;
	t149 = t133 * t160;
	t147 = t155 * t135;
	t143 = t132 * t170 + t150 * t135;
	t110 = -t123 * t174 + t183;
	t104 = t117 * t196;
	t89 = (t195 * t132 * t113 - t114 * t149) * t135;
	t88 = -t113 * t174 + t184 + (t113 * t134 - t114 * t175) * t104;
	t87 = t196 * t198 + (qJD(1) * t147 + 0.2e1 * t133 * t144) * t117;
	t1 = [t154 * t178 + (qJD(2) * t147 - t172 * t179) * t117, t87, 0, 0; (-t99 * t162 + (0.2e1 * t164 + (qJD(1) * t89 + t86) * t192) * t132) * t133 + ((-t89 * t162 + (t89 * t167 + ((0.2e1 * t132 * t188 - t92 * t149 - t195 * t169) * t113 + (t154 * t180 + t132 * t92 + (t125 * t157 - (t92 - 0.2e1 * t170) * t132) * t117) * t114) * t93 * t135) * t132) * t100 + (t89 * t165 + (-t99 + ((-t127 + t131) * t114 * t160 + t195 * t161) * t100) * qJD(1)) * t132 * t93) * t135, (-t99 * t93 * t172 + (-0.2e1 * t164 + (-qJD(2) * t88 - t86) * t192) * t135) * t134 + (t88 * t135 * t100 * t167 + ((-qJD(2) * t99 + t88 * t165) * t135 + (t88 * t172 + (-(-t104 * t171 - t133 * t87) * t114 - ((t104 * t133 - 0.1e1) * t92 + (-t104 + t133) * qJD(2)) * t113) * t178) * t100) * t93 - ((t87 - t171) * t113 + (t152 * t104 + t153) * t114) * t173 * t192) * t132, 0, 0; 0.2e1 * (t106 * t145 + t110 * t186) * t193 + (0.2e1 * t110 * t163 - t151 * t106 * t176 + t143 * t187 + (-t151 * t111 * t177 + t110 * t90 - t143 * t185 + t145 * t91) * t107) * t96, t134 * t168 * t197 + (-t172 * t197 + (t146 * t166 + ((-qJD(4) * t106 - 0.2e1 * t163) * t123 + (-t123 * t90 + (-qJD(4) * t111 + t91) * t122) * t107) * t96) * t135) * t132, 0, t166 + 0.2e1 * (-t107 * t90 * t96 + (-t107 * t193 - t96 * t190) * t111) * t111;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end