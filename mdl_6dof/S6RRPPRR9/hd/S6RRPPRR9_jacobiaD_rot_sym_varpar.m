% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR9
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
%   Wie in S6RRPPRR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:50
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR9_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR9_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (215->39), mult. (853->106), div. (126->12), fcn. (1047->9), ass. (0->54)
	t99 = sin(pkin(6));
	t93 = t99 ^ 2;
	t100 = cos(pkin(6));
	t95 = 0.1e1 / t100 ^ 2;
	t104 = cos(qJ(1));
	t98 = t104 ^ 2;
	t89 = t98 * t93 * t95 + 0.1e1;
	t102 = sin(qJ(1));
	t97 = t102 ^ 2;
	t126 = 0.1e1 / t89 ^ 2 * t97;
	t131 = t126 * t95;
	t122 = t104 * t99;
	t88 = atan2(t122, t100);
	t84 = sin(t88);
	t85 = cos(t88);
	t72 = t85 * t100 + t84 * t122;
	t67 = 0.1e1 / t72;
	t103 = cos(qJ(2));
	t118 = t104 * t103;
	t101 = sin(qJ(2));
	t121 = t102 * t101;
	t113 = t100 * t121 - t118;
	t77 = 0.1e1 / t113;
	t94 = 0.1e1 / t100;
	t68 = 0.1e1 / t72 ^ 2;
	t78 = 0.1e1 / t113 ^ 2;
	t119 = t104 * t101;
	t120 = t102 * t103;
	t81 = -t100 * t119 - t120;
	t82 = t100 * t120 + t119;
	t71 = t81 * qJD(1) - t82 * qJD(2);
	t128 = t71 * t77 * t78;
	t115 = t100 * t118;
	t70 = -qJD(1) * t115 - qJD(2) * t118 + (qJD(2) * t100 + qJD(1)) * t121;
	t129 = t70 * t78;
	t76 = t82 ^ 2;
	t75 = t76 * t78 + 0.1e1;
	t130 = (t76 * t128 - t82 * t129) / t75 ^ 2;
	t127 = t81 * t82;
	t125 = t93 * t94;
	t124 = t102 * t68;
	t123 = t104 * t68;
	t117 = qJD(1) * t104;
	t86 = 0.1e1 / t89;
	t116 = (t86 - 0.1e1) * t99;
	t114 = -0.2e1 * t94 * t131;
	t80 = t115 - t121;
	t63 = (-t104 * t85 * t86 * t125 + t84 * t116) * t102;
	t92 = t99 * t93;
	t73 = 0.1e1 / t75;
	t69 = t67 * t68;
	t66 = t97 * t93 * t68 + 0.1e1;
	t62 = qJD(1) * t63;
	t1 = [(-t86 * t94 * t99 + t92 * t114) * t117, 0, 0, 0, 0, 0; (0.2e1 * (-t104 * t67 + t63 * t124) / t66 ^ 2 * (-t62 * t69 * t97 + t117 * t124) * t93 + ((0.2e1 * t102 * t63 * t69 - t123) * t62 + (-t63 * t123 + (-t67 + (-t92 * t131 - t116) * t84 * t123 - (t93 ^ 2 * t98 * t114 + (-t126 + (0.2e1 * t97 - t98) * t86) * t125) * t68 * t85) * t102) * qJD(1)) / t66) * t99, 0, 0, 0, 0, 0; 0.2e1 * (t78 * t127 + t77 * t80) * t130 + (-(-t82 * qJD(1) + t81 * qJD(2)) * t77 - 0.2e1 * t127 * t128 + (-t80 * t71 - (t113 * qJD(1) - t80 * qJD(2)) * t82 + t81 * t70) * t78) * t73, -0.2e1 * t130 + 0.2e1 * (-t73 * t129 + (t73 * t128 - t78 * t130) * t82) * t82, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:18
	% DurationCPUTime: 0.74s
	% Computational Cost: add. (1114->72), mult. (3196->173), div. (656->14), fcn. (4222->9), ass. (0->74)
	t125 = cos(qJ(2));
	t126 = cos(qJ(1));
	t156 = cos(pkin(6));
	t141 = t126 * t156;
	t139 = t125 * t141;
	t123 = sin(qJ(2));
	t124 = sin(qJ(1));
	t152 = t124 * t123;
	t104 = -t139 + t152;
	t122 = sin(pkin(6));
	t114 = 0.1e1 / t122;
	t119 = 0.1e1 / t125;
	t144 = t104 * t114 * t119;
	t153 = t122 * t125;
	t94 = atan2(-t104, -t153);
	t92 = sin(t94);
	t93 = cos(t94);
	t101 = t104 ^ 2;
	t115 = 0.1e1 / t122 ^ 2;
	t120 = 0.1e1 / t125 ^ 2;
	t99 = t101 * t115 * t120 + 0.1e1;
	t95 = 0.1e1 / t99;
	t166 = (t144 * t93 - t92) * t95 + t92;
	t87 = -t104 * t92 - t153 * t93;
	t84 = 0.1e1 / t87;
	t116 = 0.1e1 / t124;
	t117 = 0.1e1 / t124 ^ 2;
	t85 = 0.1e1 / t87 ^ 2;
	t149 = qJD(2) * t123;
	t158 = t125 * t92;
	t163 = t104 * t93;
	t143 = t120 * t149;
	t160 = t114 * t95;
	t134 = -t123 * t141 - t124 * t125;
	t142 = t124 * t156;
	t135 = -t123 * t126 - t125 * t142;
	t90 = -qJD(1) * t135 - qJD(2) * t134;
	t77 = (t104 * t143 + t119 * t90) * t160;
	t74 = -t77 * t163 - t92 * t90 + (t149 * t93 + t158 * t77) * t122;
	t165 = t74 * t84 * t85;
	t154 = t120 * t123;
	t136 = t104 * t154 - t119 * t134;
	t78 = t136 * t160;
	t164 = t77 * t78;
	t162 = t135 * t85;
	t161 = t135 * t93;
	t159 = t119 * t95;
	t157 = t92 * t135;
	t155 = t117 * t126;
	t151 = t126 * t125;
	t150 = qJD(1) * t126;
	t102 = t135 ^ 2;
	t81 = t102 * t85 + 0.1e1;
	t138 = qJD(2) * t156 + qJD(1);
	t88 = -qJD(1) * t139 - qJD(2) * t151 + t138 * t152;
	t148 = 0.2e1 * (-t102 * t165 + t162 * t88) / t81 ^ 2;
	t147 = 0.2e1 * t165;
	t121 = t119 * t120;
	t146 = -0.2e1 * (t101 * t121 * t149 + t104 * t120 * t90) * t115 / t99 ^ 2;
	t140 = t123 * t142;
	t108 = -t140 + t151;
	t103 = t108 ^ 2;
	t100 = t103 * t115 * t117 + 0.1e1;
	t118 = t116 * t117;
	t89 = qJD(1) * t134 + qJD(2) * t135;
	t145 = 0.2e1 * (-t103 * t118 * t150 + t108 * t117 * t89) * t115 / t100 ^ 2;
	t133 = t119 * t146 + t143 * t95;
	t97 = 0.1e1 / t100;
	t91 = -qJD(1) * t140 - t124 * t149 + t138 * t151;
	t79 = 0.1e1 / t81;
	t76 = t166 * t135;
	t75 = -t78 * t163 + t92 * t134 + (t123 * t93 + t158 * t78) * t122;
	t73 = (t136 * t146 + (t90 * t154 + t119 * t91 + (-t134 * t154 + (0.2e1 * t121 * t123 ^ 2 + t119) * t104) * qJD(2)) * t95) * t114;
	t1 = [(-t133 * t135 - t159 * t88) * t114, t73, 0, 0, 0, 0; t104 * t84 * t148 + (-t90 * t84 + (t104 * t74 + t76 * t88) * t85) * t79 - (t76 * t147 * t79 + (t76 * t148 + ((t144 * t77 * t95 + t146) * t157 + ((t95 - 0.1e1) * t77 + (-t104 * t133 - t159 * t90) * t114) * t161 - t166 * t88) * t79) * t85) * t135, (-t108 * t84 - t162 * t75) * t148 + (-t75 * t135 * t147 + t89 * t84 + (-t108 * t74 + t75 * t88 + (t104 * t164 - t91) * t157 + (-t104 * t73 + t134 * t77 - t78 * t90) * t161) * t85 + ((-qJD(2) * t78 - t77) * t92 * t123 + (t73 * t92 + (qJD(2) + t164) * t93) * t125) * t122 * t162) * t79, 0, 0, 0, 0; ((t108 * t155 - t116 * t134) * t145 + (-t89 * t155 - t116 * t91 + (-t134 * t155 + (0.2e1 * t118 * t126 ^ 2 + t116) * t108) * qJD(1)) * t97) * t114, (t116 * t88 * t97 - (t117 * t150 * t97 + t116 * t145) * t135) * t114, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:18
	% DurationCPUTime: 0.74s
	% Computational Cost: add. (969->70), mult. (3196->173), div. (656->14), fcn. (4222->9), ass. (0->73)
	t120 = sin(qJ(2));
	t121 = sin(qJ(1));
	t122 = cos(qJ(2));
	t123 = cos(qJ(1));
	t151 = cos(pkin(6));
	t136 = t123 * t151;
	t103 = t120 * t136 + t121 * t122;
	t119 = sin(pkin(6));
	t111 = 0.1e1 / t119;
	t113 = 0.1e1 / t120;
	t139 = t103 * t111 * t113;
	t148 = t119 * t120;
	t96 = atan2(-t103, t148);
	t90 = sin(t96);
	t91 = cos(t96);
	t112 = 0.1e1 / t119 ^ 2;
	t114 = 0.1e1 / t120 ^ 2;
	t99 = t103 ^ 2;
	t97 = t99 * t112 * t114 + 0.1e1;
	t92 = 0.1e1 / t97;
	t162 = (t91 * t139 + t90) * t92 - t90;
	t85 = -t90 * t103 + t91 * t148;
	t82 = 0.1e1 / t85;
	t116 = 0.1e1 / t121;
	t117 = 0.1e1 / t121 ^ 2;
	t83 = 0.1e1 / t85 ^ 2;
	t144 = qJD(2) * t122;
	t153 = t120 * t90;
	t158 = t103 * t91;
	t138 = t114 * t144;
	t155 = t111 * t92;
	t137 = t121 * t151;
	t135 = t120 * t137;
	t146 = t123 * t122;
	t147 = t121 * t120;
	t89 = -qJD(1) * t135 - qJD(2) * t147 + (qJD(2) * t151 + qJD(1)) * t146;
	t75 = (t103 * t138 - t113 * t89) * t155;
	t72 = -t75 * t158 - t90 * t89 + (t91 * t144 - t75 * t153) * t119;
	t161 = t72 * t82 * t83;
	t102 = -t122 * t136 + t147;
	t150 = t114 * t122;
	t133 = t102 * t113 + t103 * t150;
	t76 = t133 * t155;
	t160 = t75 * t76;
	t115 = t113 * t114;
	t159 = (t103 * t114 * t89 - t115 * t99 * t144) * t112 / t97 ^ 2;
	t132 = t135 - t146;
	t157 = t132 * t83;
	t156 = t132 * t91;
	t154 = t113 * t92;
	t152 = t90 * t132;
	t149 = t117 * t123;
	t145 = qJD(1) * t123;
	t101 = t132 ^ 2;
	t79 = t101 * t83 + 0.1e1;
	t131 = t123 * t120 + t122 * t137;
	t87 = t103 * qJD(1) + t131 * qJD(2);
	t143 = 0.2e1 * (-t101 * t161 + t87 * t157) / t79 ^ 2;
	t142 = 0.2e1 * t161;
	t141 = -0.2e1 * t159;
	t100 = t131 ^ 2;
	t118 = t116 * t117;
	t86 = t102 * qJD(1) + t132 * qJD(2);
	t98 = t100 * t117 * t112 + 0.1e1;
	t140 = 0.2e1 * (-t100 * t118 * t145 - t117 * t131 * t86) * t112 / t98 ^ 2;
	t130 = 0.2e1 * t113 * t159 + t92 * t138;
	t94 = 0.1e1 / t98;
	t88 = t131 * qJD(1) + t103 * qJD(2);
	t77 = 0.1e1 / t79;
	t74 = t162 * t132;
	t73 = -t76 * t158 + t90 * t102 + (t122 * t91 - t76 * t153) * t119;
	t71 = (t133 * t141 + (t89 * t150 + t113 * t88 + (-t102 * t150 + (-0.2e1 * t115 * t122 ^ 2 - t113) * t103) * qJD(2)) * t92) * t111;
	t1 = [(-t130 * t132 + t87 * t154) * t111, t71, 0, 0, 0, 0; t103 * t82 * t143 + (-t89 * t82 + (t103 * t72 - t74 * t87) * t83) * t77 - (-t74 * t142 * t77 + (-t74 * t143 + ((-t75 * t92 * t139 + t141) * t152 + ((t92 - 0.1e1) * t75 + (-t130 * t103 + t89 * t154) * t111) * t156 + t162 * t87) * t77) * t83) * t132, (t131 * t82 - t73 * t157) * t143 + (-t73 * t132 * t142 + t86 * t82 + (t131 * t72 + t73 * t87 + (t103 * t160 + t88) * t152 + (t102 * t75 - t103 * t71 - t76 * t89) * t156) * t83 + ((-qJD(2) * t76 - t75) * t90 * t122 + (-t71 * t90 + (-qJD(2) - t160) * t91) * t120) * t119 * t157) * t77, 0, 0, 0, 0; ((-t102 * t116 - t131 * t149) * t140 + (-t86 * t149 + t116 * t88 + (-t102 * t149 - (0.2e1 * t118 * t123 ^ 2 + t116) * t131) * qJD(1)) * t94) * t111, (t116 * t87 * t94 - (t117 * t94 * t145 + t116 * t140) * t132) * t111, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:18
	% DurationCPUTime: 0.84s
	% Computational Cost: add. (1334->89), mult. (4303->202), div. (668->14), fcn. (5516->11), ass. (0->93)
	t175 = sin(qJ(2));
	t176 = sin(qJ(1));
	t178 = cos(qJ(2));
	t179 = cos(qJ(1));
	t226 = cos(pkin(6));
	t195 = t179 * t226;
	t158 = t175 * t195 + t176 * t178;
	t196 = t176 * t226;
	t159 = t179 * t175 + t178 * t196;
	t138 = t159 * qJD(1) + t158 * qJD(2);
	t192 = t178 * t195;
	t211 = t176 * t175;
	t157 = -t192 + t211;
	t155 = t157 ^ 2;
	t173 = sin(pkin(6));
	t169 = 0.1e1 / t173 ^ 2;
	t171 = 0.1e1 / t178 ^ 2;
	t153 = t155 * t169 * t171 + 0.1e1;
	t170 = 0.1e1 / t178;
	t172 = t170 * t171;
	t208 = qJD(2) * t175;
	t221 = (t138 * t157 * t171 + t155 * t172 * t208) * t169 / t153 ^ 2;
	t229 = -0.2e1 * t221;
	t168 = 0.1e1 / t173;
	t228 = t157 * t168;
	t213 = t173 * t178;
	t152 = atan2(t157, t213);
	t148 = sin(t152);
	t149 = cos(t152);
	t150 = 0.1e1 / t153;
	t200 = t170 * t228;
	t227 = (t149 * t200 - t148) * t150 + t148;
	t132 = t148 * t157 + t149 * t213;
	t129 = 0.1e1 / t132;
	t193 = t175 * t196;
	t210 = t179 * t178;
	t161 = -t193 + t210;
	t174 = sin(qJ(5));
	t177 = cos(qJ(5));
	t214 = t173 * t176;
	t145 = t161 * t174 + t177 * t214;
	t141 = 0.1e1 / t145;
	t130 = 0.1e1 / t132 ^ 2;
	t142 = 0.1e1 / t145 ^ 2;
	t156 = t159 ^ 2;
	t125 = t156 * t130 + 0.1e1;
	t191 = qJD(2) * t226 + qJD(1);
	t207 = qJD(2) * t178;
	t136 = -qJD(1) * t192 - t179 * t207 + t191 * t211;
	t219 = t136 * t130;
	t197 = t171 * t208;
	t186 = (t138 * t170 + t157 * t197) * t168;
	t121 = t150 * t186;
	t188 = -t148 * t213 + t149 * t157;
	t201 = t149 * t173 * t175;
	t117 = -qJD(2) * t201 + t188 * t121 + t148 * t138;
	t224 = t117 * t129 * t130;
	t225 = (-t156 * t224 - t159 * t219) / t125 ^ 2;
	t215 = t171 * t175;
	t187 = t157 * t215 + t158 * t170;
	t122 = t187 * t168 * t150;
	t118 = t188 * t122 + t148 * t158 - t201;
	t223 = t118 * t159;
	t137 = t158 * qJD(1) + t159 * qJD(2);
	t209 = qJD(1) * t173;
	t198 = t179 * t209;
	t127 = t145 * qJD(5) + t137 * t177 + t174 * t198;
	t144 = -t161 * t177 + t174 * t214;
	t140 = t144 ^ 2;
	t135 = t140 * t142 + 0.1e1;
	t218 = t142 * t144;
	t206 = qJD(5) * t144;
	t128 = -t137 * t174 + t177 * t198 - t206;
	t220 = t128 * t141 * t142;
	t222 = (t127 * t218 - t140 * t220) / t135 ^ 2;
	t217 = t148 * t159;
	t216 = t149 * t159;
	t212 = t173 * t179;
	t205 = -0.2e1 * t225;
	t204 = -0.2e1 * t224;
	t203 = 0.2e1 * t222;
	t202 = t144 * t220;
	t199 = t176 * t209;
	t194 = t170 * t229;
	t189 = -t177 * t141 - t174 * t218;
	t147 = -t158 * t174 + t177 * t212;
	t146 = t158 * t177 + t174 * t212;
	t139 = -qJD(1) * t193 - t176 * t208 + t191 * t210;
	t133 = 0.1e1 / t135;
	t123 = 0.1e1 / t125;
	t120 = t227 * t159;
	t116 = (t187 * t229 + (t138 * t215 + t139 * t170 + (t158 * t215 + (0.2e1 * t172 * t175 ^ 2 + t170) * t157) * qJD(2)) * t150) * t168;
	t1 = [(t159 * t194 + (-t136 * t170 + t159 * t197) * t150) * t168, t116, 0, 0, 0, 0; t157 * t129 * t205 + (t138 * t129 + (-t117 * t157 - t120 * t136) * t130) * t123 + ((t120 * t204 - t227 * t219) * t123 + (t120 * t205 + ((-t121 * t150 * t200 + 0.2e1 * t221) * t217 + (t194 * t228 + t121 + (-t121 + t186) * t150) * t216) * t123) * t130) * t159, 0.2e1 * (t129 * t161 - t130 * t223) * t225 + (t204 * t223 + t137 * t129 + (t161 * t117 - t118 * t136 + (-t173 * t207 + t116 * t157 + t122 * t138 + (-t122 * t213 + t158) * t121) * t216 + (-t121 * t122 * t157 + t139 + (-t116 * t178 + (qJD(2) * t122 + t121) * t175) * t173) * t217) * t130) * t123, 0, 0, 0, 0; (-t141 * t146 + t147 * t218) * t203 + ((t147 * qJD(5) + t139 * t177 - t174 * t199) * t141 + 0.2e1 * t147 * t202 + (-t146 * t128 - (-t146 * qJD(5) - t139 * t174 - t177 * t199) * t144 - t147 * t127) * t142) * t133, t189 * t159 * t203 + (t189 * t136 + ((-qJD(5) * t141 - 0.2e1 * t202) * t174 + (t127 * t174 + (-t128 + t206) * t177) * t142) * t159) * t133, 0, 0, -0.2e1 * t222 + 0.2e1 * (t127 * t142 * t133 + (-t133 * t220 - t142 * t222) * t144) * t144, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:20
	% DurationCPUTime: 2.68s
	% Computational Cost: add. (4522->145), mult. (13478->290), div. (726->12), fcn. (17045->13), ass. (0->124)
	t260 = sin(qJ(1));
	t339 = cos(pkin(6));
	t340 = sin(qJ(2));
	t288 = t339 * t340;
	t252 = t260 * t288;
	t263 = cos(qJ(2));
	t341 = cos(qJ(1));
	t289 = t339 * t341;
	t298 = t341 * qJD(1);
	t299 = qJD(2) * t340;
	t257 = sin(pkin(6));
	t303 = t257 * t341;
	t351 = (qJD(2) * t289 + t298) * t263 - qJD(1) * t252 - t260 * t299 + qJD(5) * t303;
	t282 = t341 * t288;
	t318 = t260 * t263;
	t246 = t282 + t318;
	t259 = sin(qJ(5));
	t262 = cos(qJ(5));
	t280 = t246 * t262 + t259 * t303;
	t228 = t280 ^ 2;
	t302 = t257 * t340;
	t276 = -t339 * t259 + t262 * t302;
	t241 = 0.1e1 / t276 ^ 2;
	t219 = t228 * t241 + 0.1e1;
	t217 = 0.1e1 / t219;
	t244 = t259 * t302 + t339 * t262;
	t319 = t257 * t263;
	t300 = qJD(2) * t319;
	t230 = t244 * qJD(5) - t262 * t300;
	t240 = 0.1e1 / t276;
	t325 = t280 * t241;
	t316 = qJD(1) * t260;
	t342 = t246 * qJD(5) + t257 * t316;
	t349 = -t342 * t259 + t351 * t262;
	t285 = -t230 * t325 - t240 * t349;
	t186 = t285 * t217;
	t220 = atan2(t280, -t276);
	t215 = sin(t220);
	t216 = cos(t220);
	t286 = t215 * t276 + t216 * t280;
	t181 = t286 * t186 + t215 * t349 + t216 * t230;
	t198 = t215 * t280 - t216 * t276;
	t196 = 0.1e1 / t198 ^ 2;
	t350 = t181 * t196;
	t203 = t351 * t259 + t342 * t262;
	t287 = t341 * t263 - t252;
	t320 = t257 * t260;
	t232 = t287 * t259 + t262 * t320;
	t258 = sin(qJ(6));
	t261 = cos(qJ(6));
	t275 = t339 * t318 + t341 * t340;
	t322 = t275 * t261;
	t211 = t232 * t258 + t322;
	t348 = 0.2e1 * t211;
	t195 = 0.1e1 / t198;
	t347 = t195 * t350;
	t231 = t259 * t320 - t287 * t262;
	t297 = 0.2e1 * t231 * t347;
	t292 = t257 * t298;
	t294 = -qJD(1) * t282 - t275 * qJD(2) - t263 * t316;
	t205 = t232 * qJD(5) + t259 * t292 - t294 * t262;
	t331 = t205 * t196;
	t346 = -t331 + t297;
	t344 = t230 * t241;
	t245 = t260 * t340 - t263 * t289;
	t304 = t280 * t319;
	t281 = t240 * t245 + t241 * t304;
	t343 = t262 * t281;
	t212 = t232 * t261 - t258 * t275;
	t208 = 0.1e1 / t212;
	t209 = 0.1e1 / t212 ^ 2;
	t227 = t231 ^ 2;
	t194 = t227 * t196 + 0.1e1;
	t338 = (-t227 * t347 + t231 * t331) / t194 ^ 2;
	t206 = -t231 * qJD(5) + t294 * t259 + t262 * t292;
	t225 = t245 * qJD(1) - t287 * qJD(2);
	t190 = t212 * qJD(6) + t206 * t258 - t225 * t261;
	t207 = t211 ^ 2;
	t201 = t207 * t209 + 0.1e1;
	t330 = t209 * t211;
	t314 = qJD(6) * t211;
	t191 = t206 * t261 + t225 * t258 - t314;
	t334 = t191 * t208 * t209;
	t337 = (t190 * t330 - t207 * t334) / t201 ^ 2;
	t327 = t240 * t344;
	t335 = (t228 * t327 + t325 * t349) / t219 ^ 2;
	t333 = t196 * t231;
	t199 = 0.1e1 / t201;
	t332 = t199 * t209;
	t329 = t215 * t231;
	t328 = t216 * t231;
	t326 = t280 * t240;
	t324 = t280 * t244;
	t323 = t275 * t259;
	t321 = t275 * t262;
	t315 = qJD(5) * t259;
	t312 = 0.2e1 * t338;
	t311 = -0.2e1 * t337;
	t310 = -0.2e1 * t335;
	t309 = 0.2e1 * t335;
	t307 = t209 * t337;
	t306 = t190 * t332;
	t305 = t211 * t334;
	t296 = t240 * t309;
	t295 = 0.2e1 * t305;
	t279 = -t246 * t259 + t262 * t303;
	t214 = t245 * t258 + t261 * t279;
	t213 = -t245 * t261 + t258 * t279;
	t284 = -t258 * t208 + t261 * t330;
	t283 = t240 * t279 + t241 * t324;
	t278 = -t215 + (t216 * t326 + t215) * t217;
	t277 = -qJD(6) * t323 + t294;
	t273 = -qJD(5) * t321 - t287 * qJD(6) + t225 * t259;
	t229 = t276 * qJD(5) + t259 * t300;
	t226 = t275 * qJD(1) + t246 * qJD(2);
	t222 = -t287 * t258 - t259 * t322;
	t192 = 0.1e1 / t194;
	t189 = t217 * t343;
	t188 = t283 * t217;
	t183 = (-t215 * t245 - t216 * t319) * t262 + t286 * t189;
	t182 = -t286 * t188 + t215 * t279 + t216 * t244;
	t180 = t283 * t309 + (-0.2e1 * t324 * t327 + t203 * t240 + (-t229 * t280 - t230 * t279 - t244 * t349) * t241) * t217;
	t178 = t310 * t343 + (-t281 * t315 + (0.2e1 * t304 * t327 + t226 * t240 + (t230 * t245 + (t263 * t349 - t280 * t299) * t257) * t241) * t262) * t217;
	t1 = [-t231 * t296 + (t205 * t240 + t231 * t344) * t217, t178, 0, 0, t180, 0; -0.2e1 * t280 * t195 * t338 + (t349 * t195 - t280 * t350 - (t278 * t205 + ((-t186 * t217 * t326 + t310) * t215 + (-t280 * t296 - t186 + (t186 - t285) * t217) * t216) * t231) * t333) * t192 + (t346 * t192 + t333 * t312) * t278 * t231, (t183 * t333 - t195 * t321) * t312 + ((-t225 * t262 - t275 * t315) * t195 + t346 * t183 + (-t321 * t181 - (t178 * t280 + t189 * t349 + (t262 * t299 + t263 * t315) * t257 + (t189 * t276 - t245 * t262) * t186) * t328 - (t245 * t315 + t178 * t276 - t189 * t230 - t226 * t262 + (-t189 * t280 + t262 * t319) * t186) * t329) * t196) * t192, 0, 0, (t182 * t333 - t195 * t232) * t312 + (t182 * t297 + t206 * t195 + (-t232 * t181 - t182 * t205 - (t180 * t280 - t188 * t349 + t229 + (-t188 * t276 + t279) * t186) * t328 - (t180 * t276 + t188 * t230 - t203 + (t188 * t280 - t244) * t186) * t329) * t196) * t192, 0; 0.2e1 * (-t208 * t213 + t214 * t330) * t337 + ((qJD(6) * t214 - t203 * t258 - t226 * t261) * t208 + t214 * t295 + (-t213 * t191 - (-qJD(6) * t213 - t203 * t261 + t226 * t258) * t211 - t214 * t190) * t209) * t199, (t307 * t348 - t306) * t222 + (-t191 * t332 + t208 * t311) * (-t258 * t323 + t287 * t261) + ((t258 * t273 + t261 * t277) * t208 - (-t258 * t277 + t261 * t273) * t330 + t222 * t295) * t199, 0, 0, t284 * t231 * t311 + (t284 * t205 + ((-qJD(6) * t208 - 0.2e1 * t305) * t261 + (t190 * t261 + (t191 - t314) * t258) * t209) * t231) * t199, t311 + (t306 + (-t199 * t334 - t307) * t211) * t348;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end