% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR11
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
%   Wie in S6RRPPRR11_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:53
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR11_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR11_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.41s
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
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.77s
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
	t166 = (t93 * t144 - t92) * t95 + t92;
	t87 = -t92 * t104 - t93 * t153;
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
	t135 = -t126 * t123 - t125 * t142;
	t90 = -t135 * qJD(1) - t134 * qJD(2);
	t77 = (t104 * t143 + t119 * t90) * t160;
	t74 = -t77 * t163 - t92 * t90 + (t93 * t149 + t77 * t158) * t122;
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
	t148 = 0.2e1 * (-t102 * t165 + t88 * t162) / t81 ^ 2;
	t147 = 0.2e1 * t165;
	t121 = t119 * t120;
	t146 = -0.2e1 * (t101 * t121 * t149 + t104 * t120 * t90) * t115 / t99 ^ 2;
	t140 = t123 * t142;
	t108 = -t140 + t151;
	t103 = t108 ^ 2;
	t100 = t103 * t117 * t115 + 0.1e1;
	t118 = t116 * t117;
	t89 = t134 * qJD(1) + t135 * qJD(2);
	t145 = 0.2e1 * (-t103 * t118 * t150 + t108 * t117 * t89) * t115 / t100 ^ 2;
	t133 = t119 * t146 + t95 * t143;
	t97 = 0.1e1 / t100;
	t91 = -qJD(1) * t140 - t124 * t149 + t138 * t151;
	t79 = 0.1e1 / t81;
	t76 = t166 * t135;
	t75 = -t78 * t163 + t92 * t134 + (t123 * t93 + t78 * t158) * t122;
	t73 = (t136 * t146 + (t90 * t154 + t119 * t91 + (-t134 * t154 + (0.2e1 * t121 * t123 ^ 2 + t119) * t104) * qJD(2)) * t95) * t114;
	t1 = [(-t133 * t135 - t88 * t159) * t114, t73, 0, 0, 0, 0; t104 * t84 * t148 + (-t90 * t84 + (t104 * t74 + t76 * t88) * t85) * t79 - (t76 * t147 * t79 + (t76 * t148 + ((t77 * t95 * t144 + t146) * t157 + ((t95 - 0.1e1) * t77 + (-t133 * t104 - t90 * t159) * t114) * t161 - t166 * t88) * t79) * t85) * t135, (-t108 * t84 - t75 * t162) * t148 + (-t75 * t135 * t147 + t89 * t84 + (-t108 * t74 + t75 * t88 + (t104 * t164 - t91) * t157 + (-t104 * t73 + t134 * t77 - t78 * t90) * t161) * t85 + ((-qJD(2) * t78 - t77) * t92 * t123 + (t73 * t92 + (qJD(2) + t164) * t93) * t125) * t122 * t162) * t79, 0, 0, 0, 0; ((t108 * t155 - t116 * t134) * t145 + (-t89 * t155 - t116 * t91 + (-t134 * t155 + (0.2e1 * t118 * t126 ^ 2 + t116) * t108) * qJD(1)) * t97) * t114, (t116 * t88 * t97 - (t117 * t97 * t150 + t116 * t145) * t135) * t114, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.77s
	% Computational Cost: add. (1150->84), mult. (3755->189), div. (644->14), fcn. (4884->11), ass. (0->87)
	t143 = sin(qJ(2));
	t144 = cos(qJ(2));
	t145 = cos(qJ(1));
	t181 = cos(pkin(6));
	t160 = t145 * t181;
	t188 = sin(qJ(1));
	t125 = t143 * t160 + t188 * t144;
	t156 = t181 * t188;
	t127 = t145 * t143 + t144 * t156;
	t105 = t125 * qJD(1) + t127 * qJD(2);
	t155 = t143 * t156;
	t173 = t144 * t145;
	t128 = -t155 + t173;
	t123 = t128 ^ 2;
	t141 = sin(pkin(6));
	t176 = t141 * t143;
	t120 = atan2(-t125, t176);
	t116 = sin(t120);
	t117 = cos(t120);
	t98 = -t116 * t125 + t117 * t176;
	t96 = 0.1e1 / t98 ^ 2;
	t182 = t128 * t96;
	t164 = t188 * t143;
	t107 = -qJD(1) * t155 - qJD(2) * t164 + (qJD(2) * t181 + qJD(1)) * t173;
	t172 = qJD(2) * t144;
	t137 = 0.1e1 / t143;
	t138 = 0.1e1 / t143 ^ 2;
	t162 = t138 * t172;
	t152 = -t107 * t137 + t125 * t162;
	t122 = t125 ^ 2;
	t136 = 0.1e1 / t141 ^ 2;
	t121 = t122 * t136 * t138 + 0.1e1;
	t118 = 0.1e1 / t121;
	t135 = 0.1e1 / t141;
	t178 = t118 * t135;
	t89 = t152 * t178;
	t85 = (-t125 * t89 + t141 * t172) * t117 + (-t89 * t176 - t107) * t116;
	t95 = 0.1e1 / t98;
	t97 = t95 * t96;
	t186 = t85 * t97;
	t93 = t123 * t96 + 0.1e1;
	t171 = 0.2e1 * (-t105 * t182 - t123 * t186) / t93 ^ 2;
	t166 = t125 * t135 * t137;
	t190 = (t117 * t166 + t116) * t118 - t116;
	t140 = sin(pkin(11));
	t142 = cos(pkin(11));
	t165 = t141 * t188;
	t113 = t127 * t140 + t142 * t165;
	t109 = 0.1e1 / t113;
	t110 = 0.1e1 / t113 ^ 2;
	t189 = 0.2e1 * t128;
	t185 = t105 * t96;
	t139 = t137 * t138;
	t184 = 0.1e1 / t121 ^ 2 * (t107 * t125 * t138 - t122 * t139 * t172) * t136;
	t157 = t144 * t160;
	t124 = t164 - t157;
	t177 = t138 * t144;
	t153 = t124 * t137 + t125 * t177;
	t90 = t153 * t178;
	t183 = t125 * t90;
	t161 = t188 * qJD(1);
	t104 = -qJD(1) * t157 - t145 * t172 + (qJD(2) * t156 + t161) * t143;
	t175 = t141 * t145;
	t163 = qJD(1) * t175;
	t103 = -t104 * t140 + t142 * t163;
	t180 = t103 * t109 * t110;
	t112 = -t127 * t142 + t140 * t165;
	t179 = t110 * t112;
	t174 = t142 * t109;
	t108 = t112 ^ 2;
	t101 = t108 * t110 + 0.1e1;
	t102 = t104 * t142 + t140 * t163;
	t170 = 0.2e1 / t101 ^ 2 * (t102 * t179 - t108 * t180);
	t169 = -0.2e1 * t184;
	t168 = t116 * t182;
	t167 = t117 * t182;
	t159 = 0.2e1 * t112 * t180;
	t158 = t141 * t161;
	t115 = -t124 * t140 + t142 * t175;
	t114 = t124 * t142 + t140 * t175;
	t106 = t127 * qJD(1) + qJD(2) * t125;
	t99 = 0.1e1 / t101;
	t91 = 0.1e1 / t93;
	t88 = t190 * t128;
	t86 = (t141 * t144 - t183) * t117 + (-t90 * t176 + t124) * t116;
	t84 = (t153 * t169 + (t107 * t177 + t106 * t137 + (-t124 * t177 + (-0.2e1 * t139 * t144 ^ 2 - t137) * t125) * qJD(2)) * t118) * t135;
	t1 = [(t137 * t184 * t189 + (t105 * t137 + t128 * t162) * t118) * t135, t84, 0, 0, 0, 0; (-t107 * t95 + (t105 * t88 + t125 * t85) * t96 + (0.2e1 * t186 * t88 - (-t118 * t89 * t166 + t169) * t168 - (t166 * t169 - t89 + (-t152 * t135 + t89) * t118) * t167 + t190 * t185) * t128) * t91 + (t125 * t95 + t88 * t182) * t171, (t127 * t95 + t86 * t182) * t171 + (t86 * t185 + t104 * t95 + (t86 * t97 * t189 + t127 * t96) * t85 - (-t107 * t90 + t124 * t89 - t125 * t84 + (-t89 * t90 - qJD(2)) * t176) * t167 - (t89 * t183 + t106 + (-t143 * t84 + (-qJD(2) * t90 - t89) * t144) * t141) * t168) * t91, 0, 0, 0, 0; (-t109 * t114 + t115 * t179) * t170 + ((t106 * t142 - t140 * t158) * t109 + t115 * t159 + (-t114 * t103 - (-t106 * t140 - t142 * t158) * t112 - t115 * t102) * t110) * t99, (t140 * t179 + t174) * t99 * t105 + (t170 * t174 + t140 * t99 * t159 + (t103 * t142 * t99 + (-t102 * t99 + t112 * t170) * t140) * t110) * t128, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:52
	% DurationCPUTime: 0.92s
	% Computational Cost: add. (1643->91), mult. (4303->199), div. (668->14), fcn. (5516->11), ass. (0->92)
	t181 = sin(qJ(2));
	t182 = sin(qJ(1));
	t183 = cos(qJ(2));
	t184 = cos(qJ(1));
	t231 = cos(pkin(6));
	t200 = t184 * t231;
	t162 = t181 * t200 + t182 * t183;
	t180 = sin(pkin(6));
	t219 = t180 * t181;
	t156 = atan2(-t162, t219);
	t152 = sin(t156);
	t153 = cos(t156);
	t159 = t162 ^ 2;
	t175 = 0.1e1 / t180 ^ 2;
	t178 = 0.1e1 / t181 ^ 2;
	t157 = t159 * t175 * t178 + 0.1e1;
	t154 = 0.1e1 / t157;
	t174 = 0.1e1 / t180;
	t177 = 0.1e1 / t181;
	t205 = t162 * t174 * t177;
	t232 = t154 * (t153 * t205 + t152) - t152;
	t136 = -t152 * t162 + t153 * t219;
	t133 = 0.1e1 / t136;
	t201 = t182 * t231;
	t164 = t184 * t181 + t183 * t201;
	t176 = pkin(11) + qJ(5);
	t172 = sin(t176);
	t173 = cos(t176);
	t218 = t180 * t182;
	t149 = t164 * t172 + t173 * t218;
	t145 = 0.1e1 / t149;
	t134 = 0.1e1 / t136 ^ 2;
	t146 = 0.1e1 / t149 ^ 2;
	t196 = qJD(2) * t231 + qJD(1);
	t198 = t181 * t201;
	t213 = qJD(2) * t181;
	t215 = t184 * t183;
	t143 = -qJD(1) * t198 - t182 * t213 + t196 * t215;
	t212 = qJD(2) * t183;
	t202 = t178 * t212;
	t191 = -t143 * t177 + t162 * t202;
	t221 = t154 * t174;
	t125 = t191 * t221;
	t193 = -t152 * t219 - t153 * t162;
	t206 = t153 * t180 * t183;
	t121 = qJD(2) * t206 + t125 * t193 - t152 * t143;
	t230 = t121 * t133 * t134;
	t197 = t183 * t200;
	t216 = t182 * t181;
	t161 = -t197 + t216;
	t220 = t178 * t183;
	t192 = t161 * t177 + t162 * t220;
	t126 = t192 * t221;
	t122 = t126 * t193 + t152 * t161 + t206;
	t165 = -t198 + t215;
	t229 = t122 * t165;
	t140 = -qJD(1) * t197 - t184 * t212 + t196 * t216;
	t214 = qJD(1) * t180;
	t203 = t184 * t214;
	t130 = qJD(5) * t149 + t140 * t173 + t172 * t203;
	t148 = -t164 * t173 + t172 * t218;
	t144 = t148 ^ 2;
	t139 = t144 * t146 + 0.1e1;
	t224 = t146 * t148;
	t211 = qJD(5) * t148;
	t131 = -t140 * t172 + t173 * t203 - t211;
	t227 = t131 * t145 * t146;
	t228 = (t130 * t224 - t144 * t227) / t139 ^ 2;
	t179 = t177 * t178;
	t226 = (t143 * t162 * t178 - t159 * t179 * t212) * t175 / t157 ^ 2;
	t141 = t162 * qJD(1) + qJD(2) * t164;
	t225 = t141 * t134;
	t223 = t152 * t165;
	t222 = t153 * t165;
	t217 = t180 * t184;
	t160 = t165 ^ 2;
	t129 = t160 * t134 + 0.1e1;
	t210 = 0.2e1 * (-t160 * t230 - t165 * t225) / t129 ^ 2;
	t209 = 0.2e1 * t230;
	t208 = 0.2e1 * t228;
	t207 = -0.2e1 * t226;
	t204 = t182 * t214;
	t199 = 0.2e1 * t148 * t227;
	t194 = t173 * t145 + t172 * t224;
	t151 = -t161 * t172 + t173 * t217;
	t150 = t161 * t173 + t172 * t217;
	t142 = qJD(1) * t164 + qJD(2) * t162;
	t137 = 0.1e1 / t139;
	t127 = 0.1e1 / t129;
	t124 = t232 * t165;
	t120 = (t192 * t207 + (t143 * t220 + t142 * t177 + (-t161 * t220 + (-0.2e1 * t179 * t183 ^ 2 - t177) * t162) * qJD(2)) * t154) * t174;
	t1 = [(0.2e1 * t165 * t177 * t226 + (t141 * t177 + t165 * t202) * t154) * t174, t120, 0, 0, 0, 0; t162 * t133 * t210 + (-t143 * t133 + (t121 * t162 + t124 * t141) * t134) * t127 + ((t124 * t209 + t232 * t225) * t127 + (t124 * t210 + (-(-t125 * t154 * t205 + t207) * t223 - (t205 * t207 - t125 + (-t174 * t191 + t125) * t154) * t222) * t127) * t134) * t165, (t133 * t164 + t134 * t229) * t210 + (t209 * t229 + t140 * t133 + (t164 * t121 + t122 * t141 - (-t180 * t213 - t120 * t162 - t126 * t143 + (-t126 * t219 + t161) * t125) * t222 - (t125 * t126 * t162 + t142 + (-t120 * t181 + (-qJD(2) * t126 - t125) * t183) * t180) * t223) * t134) * t127, 0, 0, 0, 0; (-t145 * t150 + t151 * t224) * t208 + ((qJD(5) * t151 + t142 * t173 - t172 * t204) * t145 + t151 * t199 + (-t150 * t131 - (-qJD(5) * t150 - t142 * t172 - t173 * t204) * t148 - t151 * t130) * t146) * t137, t194 * t165 * t208 + (t194 * t141 + ((qJD(5) * t145 + t199) * t172 + (-t130 * t172 + (t131 - t211) * t173) * t146) * t165) * t137, 0, 0, -0.2e1 * t228 + 0.2e1 * (t130 * t146 * t137 + (-t137 * t227 - t146 * t228) * t148) * t148, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:54
	% DurationCPUTime: 2.78s
	% Computational Cost: add. (8428->148), mult. (13478->295), div. (726->12), fcn. (17045->13), ass. (0->124)
	t266 = sin(qJ(2));
	t267 = sin(qJ(1));
	t269 = cos(qJ(2));
	t270 = cos(qJ(1));
	t348 = cos(pkin(6));
	t300 = t270 * t348;
	t285 = -t266 * t300 - t267 * t269;
	t301 = t267 * t348;
	t287 = t270 * t266 + t269 * t301;
	t264 = sin(pkin(6));
	t326 = t264 * t270;
	t361 = t287 * qJD(1) - t285 * qJD(2) + qJD(5) * t326;
	t252 = t267 * t266 - t269 * t300;
	t263 = pkin(11) + qJ(5);
	t261 = sin(t263);
	t262 = cos(t263);
	t289 = t252 * t262 + t261 * t326;
	t237 = t289 ^ 2;
	t327 = t264 * t269;
	t286 = -t348 * t261 - t262 * t327;
	t247 = 0.1e1 / t286 ^ 2;
	t226 = t237 * t247 + 0.1e1;
	t220 = 0.1e1 / t226;
	t250 = -t261 * t327 + t348 * t262;
	t322 = qJD(2) * t266;
	t305 = t264 * t322;
	t235 = t250 * qJD(5) - t262 * t305;
	t246 = 0.1e1 / t286;
	t334 = t289 * t247;
	t323 = qJD(1) * t267;
	t349 = t252 * qJD(5) + t264 * t323;
	t359 = -t349 * t261 + t361 * t262;
	t293 = -t235 * t334 - t246 * t359;
	t193 = t293 * t220;
	t227 = atan2(t289, -t286);
	t214 = sin(t227);
	t215 = cos(t227);
	t295 = t214 * t286 + t215 * t289;
	t188 = t295 * t193 + t214 * t359 + t215 * t235;
	t205 = t214 * t289 - t215 * t286;
	t203 = 0.1e1 / t205 ^ 2;
	t360 = t188 * t203;
	t210 = t361 * t261 + t349 * t262;
	t328 = t264 * t267;
	t239 = t287 * t261 + t262 * t328;
	t296 = t266 * t301;
	t325 = t270 * t269;
	t254 = -t296 + t325;
	t265 = sin(qJ(6));
	t268 = cos(qJ(6));
	t222 = t239 * t265 - t254 * t268;
	t358 = 0.2e1 * t222;
	t202 = 0.1e1 / t205;
	t357 = t202 * t360;
	t350 = -t261 * t328 + t287 * t262;
	t297 = -0.2e1 * t350 * t357;
	t281 = (t348 * qJD(1) + qJD(2)) * t325 - qJD(2) * t296 - t266 * t323;
	t306 = qJD(1) * t326;
	t212 = t239 * qJD(5) + t261 * t306 - t281 * t262;
	t340 = t212 * t203;
	t356 = -t340 + t297;
	t354 = t235 * t247;
	t329 = t264 * t266;
	t309 = t289 * t329;
	t288 = t246 * t285 + t247 * t309;
	t353 = t262 * t288;
	t352 = -t254 * t261 * qJD(6) - t281;
	t332 = t254 * t262;
	t351 = qJD(5) * t332 - t287 * qJD(6);
	t223 = t239 * t268 + t254 * t265;
	t217 = 0.1e1 / t223;
	t218 = 0.1e1 / t223 ^ 2;
	t236 = t350 ^ 2;
	t199 = t236 * t203 + 0.1e1;
	t347 = (-t236 * t357 - t340 * t350) / t199 ^ 2;
	t213 = t350 * qJD(5) + t281 * t261 + t262 * t306;
	t232 = t285 * qJD(1) - t287 * qJD(2);
	t200 = t223 * qJD(6) + t213 * t265 - t232 * t268;
	t216 = t222 ^ 2;
	t208 = t216 * t218 + 0.1e1;
	t337 = t218 * t222;
	t319 = qJD(6) * t222;
	t201 = t213 * t268 + t232 * t265 - t319;
	t343 = t201 * t217 * t218;
	t346 = (t200 * t337 - t216 * t343) / t208 ^ 2;
	t336 = t246 * t354;
	t344 = (t237 * t336 + t334 * t359) / t226 ^ 2;
	t342 = t203 * t350;
	t206 = 0.1e1 / t208;
	t341 = t206 * t218;
	t339 = t214 * t350;
	t338 = t215 * t350;
	t335 = t289 * t246;
	t333 = t289 * t250;
	t331 = t261 * t265;
	t330 = t261 * t268;
	t321 = qJD(2) * t269;
	t320 = qJD(5) * t261;
	t317 = 0.2e1 * t347;
	t316 = -0.2e1 * t346;
	t315 = -0.2e1 * t344;
	t314 = 0.2e1 * t344;
	t312 = t218 * t346;
	t311 = t200 * t341;
	t310 = t222 * t343;
	t299 = t246 * t314;
	t298 = 0.2e1 * t310;
	t290 = -t252 * t261 + t262 * t326;
	t225 = t265 * t285 + t268 * t290;
	t224 = t265 * t290 - t268 * t285;
	t292 = -t265 * t217 + t268 * t337;
	t291 = t246 * t290 + t247 * t333;
	t284 = -t214 + (t215 * t335 + t214) * t220;
	t234 = t286 * qJD(5) + t261 * t305;
	t233 = -qJD(1) * t296 - t267 * t322 + (qJD(2) * t348 + qJD(1)) * t325;
	t229 = t254 * t330 - t287 * t265;
	t197 = 0.1e1 / t199;
	t196 = t220 * t353;
	t194 = t291 * t220;
	t190 = (-t214 * t285 - t215 * t329) * t262 + t295 * t196;
	t189 = -t295 * t194 + t214 * t290 + t215 * t250;
	t187 = t291 * t314 + (-0.2e1 * t333 * t336 + t210 * t246 + (-t234 * t289 - t235 * t290 - t250 * t359) * t247) * t220;
	t185 = t315 * t353 + (-t288 * t320 + (0.2e1 * t309 * t336 - t233 * t246 + (t235 * t285 + (t266 * t359 + t289 * t321) * t264) * t247) * t262) * t220;
	t1 = [t350 * t299 + (t212 * t246 - t350 * t354) * t220, t185, 0, 0, t187, 0; -0.2e1 * t289 * t202 * t347 + (t359 * t202 - t289 * t360 + (t284 * t212 - ((-t193 * t220 * t335 + t315) * t214 + (-t289 * t299 - t193 + (t193 - t293) * t220) * t215) * t350) * t342) * t197 - (t356 * t197 - t342 * t317) * t284 * t350, (-t190 * t342 + t202 * t332) * t317 + ((-t232 * t262 + t254 * t320) * t202 + t356 * t190 + (t332 * t188 + (t185 * t289 + t196 * t359 + (-t262 * t321 + t266 * t320) * t264 + (t196 * t286 - t262 * t285) * t193) * t338 + (t285 * t320 + t185 * t286 - t196 * t235 + t233 * t262 + (-t196 * t289 + t262 * t329) * t193) * t339) * t203) * t197, 0, 0, (-t189 * t342 - t202 * t239) * t317 + (t189 * t297 + t213 * t202 + (-t239 * t188 - t189 * t212 + (t187 * t289 - t194 * t359 + t234 + (-t194 * t286 + t290) * t193) * t338 + (t187 * t286 + t194 * t235 - t210 + (t194 * t289 - t250) * t193) * t339) * t203) * t197, 0; 0.2e1 * (-t217 * t224 + t225 * t337) * t346 + ((t225 * qJD(6) - t210 * t265 + t233 * t268) * t217 + t225 * t298 + (-t224 * t201 - (-t224 * qJD(6) - t210 * t268 - t233 * t265) * t222 - t225 * t200) * t218) * t206, (t312 * t358 - t311) * t229 + (-t201 * t341 + t217 * t316) * (t254 * t331 + t287 * t268) + (t229 * t298 + (t331 * t217 - t330 * t337) * t232 + (-t352 * t217 - t351 * t337) * t268 + (t351 * t217 - t352 * t337) * t265) * t206, 0, 0, -t292 * t350 * t316 + (t292 * t212 - ((-qJD(6) * t217 - 0.2e1 * t310) * t268 + (t200 * t268 + (t201 - t319) * t265) * t218) * t350) * t206, t316 + (t311 + (-t206 * t343 - t312) * t222) * t358;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end