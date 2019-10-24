% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRRP6
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
%   Wie in S5PRRRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:34
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRRRP6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:34:31
	% EndTime: 2019-10-24 10:34:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:34:31
	% EndTime: 2019-10-24 10:34:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:34:31
	% EndTime: 2019-10-24 10:34:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:34:31
	% EndTime: 2019-10-24 10:34:31
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (429->41), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->59)
	t97 = sin(qJ(2));
	t90 = t97 ^ 2;
	t99 = cos(qJ(2));
	t123 = t90 / t99 ^ 2;
	t111 = 0.1e1 + t123;
	t94 = sin(pkin(8));
	t135 = t111 * t94;
	t134 = qJD(2) * t97;
	t95 = cos(pkin(8));
	t120 = t95 * t99;
	t96 = sin(qJ(3));
	t98 = cos(qJ(3));
	t77 = t98 * t120 + t94 * t96;
	t121 = t94 * t97;
	t80 = atan2(-t121, -t99);
	t78 = sin(t80);
	t79 = cos(t80);
	t63 = -t78 * t121 - t79 * t99;
	t60 = 0.1e1 / t63;
	t73 = 0.1e1 / t77;
	t124 = t78 * t99;
	t108 = -t94 * t124 + t79 * t97;
	t115 = t79 * t121;
	t109 = -t115 + t124;
	t87 = t94 ^ 2;
	t83 = t87 * t123 + 0.1e1;
	t81 = 0.1e1 / t83;
	t67 = t81 * t135;
	t54 = t109 * t67 + t108;
	t132 = 0.2e1 * t54;
	t61 = 0.1e1 / t63 ^ 2;
	t74 = 0.1e1 / t77 ^ 2;
	t118 = qJD(2) * t99;
	t127 = t61 * t97;
	t64 = qJD(2) * t67;
	t53 = t108 * qJD(2) + t109 * t64;
	t130 = t53 * t60 * t61;
	t88 = t95 ^ 2;
	t59 = t88 * t90 * t61 + 0.1e1;
	t131 = (t118 * t127 - t90 * t130) * t88 / t59 ^ 2;
	t76 = t96 * t120 - t94 * t98;
	t125 = t74 * t76;
	t113 = t95 * t134;
	t117 = qJD(3) * t76;
	t70 = -t98 * t113 - t117;
	t126 = t70 * t73 * t74;
	t72 = t76 ^ 2;
	t68 = t72 * t74 + 0.1e1;
	t69 = -t77 * qJD(3) + t96 * t113;
	t129 = (-t69 * t125 - t72 * t126) / t68 ^ 2;
	t57 = 0.1e1 / t59;
	t128 = t57 * t61;
	t119 = t67 - t94;
	t116 = -0.2e1 * t129;
	t112 = t67 * t94 - 0.1e1;
	t110 = t98 * t125 - t73 * t96;
	t65 = 0.1e1 / t68;
	t56 = 0.2e1 * (t81 - t111 / t83 ^ 2 * t87) / t99 * t134 * t135;
	t1 = [0, t56, 0, 0, 0; 0, ((-0.2e1 * t60 * t131 + (-qJD(2) * t54 - t53) * t128) * t99 + (t61 * t131 * t132 + (t56 * t61 * t115 + t130 * t132 - qJD(2) * t60 - (-t119 * qJD(2) + t112 * t64) * t78 * t127) * t57 - (t56 * t78 + (-t112 * qJD(2) + t119 * t64) * t79) * t99 * t128) * t97) * t95, 0, 0, 0; 0, (t110 * t97 * t116 + (t110 * t118 + ((-qJD(3) * t73 - 0.2e1 * t76 * t126) * t98 + (-t69 * t98 + (t70 - t117) * t96) * t74) * t97) * t65) * t95, t116 + 0.2e1 * (-t65 * t69 * t74 + (-t65 * t126 - t74 * t129) * t76) * t76, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:34:31
	% EndTime: 2019-10-24 10:34:31
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (788->45), mult. (1346->105), div. (247->12), fcn. (1511->9), ass. (0->61)
	t125 = sin(pkin(8));
	t127 = sin(qJ(2));
	t120 = t127 ^ 2;
	t128 = cos(qJ(2));
	t149 = t120 / t128 ^ 2;
	t138 = 0.1e1 + t149;
	t163 = t125 * t138;
	t162 = qJD(2) * t127;
	t124 = qJ(3) + qJ(4);
	t114 = sin(t124);
	t115 = cos(t124);
	t126 = cos(pkin(8));
	t146 = t126 * t128;
	t104 = t125 * t114 + t115 * t146;
	t147 = t125 * t127;
	t108 = atan2(-t147, -t128);
	t106 = sin(t108);
	t107 = cos(t108);
	t95 = -t106 * t147 - t107 * t128;
	t92 = 0.1e1 / t95;
	t116 = t125 ^ 2;
	t112 = t116 * t149 + 0.1e1;
	t110 = 0.1e1 / t112;
	t97 = t110 * t163;
	t140 = t125 * t97 - 0.1e1;
	t150 = t107 * t127;
	t151 = t106 * t128;
	t153 = t125 - t97;
	t82 = -t140 * t150 - t153 * t151;
	t160 = 0.2e1 * t82;
	t100 = 0.1e1 / t104;
	t101 = 0.1e1 / t104 ^ 2;
	t93 = 0.1e1 / t95 ^ 2;
	t117 = t126 ^ 2;
	t145 = qJD(2) * t128;
	t154 = t127 * t93;
	t143 = t107 * t147;
	t96 = qJD(2) * t97;
	t81 = (-t143 + t151) * t96 + (-t125 * t151 + t150) * qJD(2);
	t157 = t81 * t92 * t93;
	t86 = t117 * t120 * t93 + 0.1e1;
	t159 = (-t120 * t157 + t145 * t154) * t117 / t86 ^ 2;
	t142 = t114 * t146;
	t103 = -t125 * t115 + t142;
	t152 = t101 * t103;
	t118 = qJD(3) + qJD(4);
	t139 = t126 * t162;
	t91 = -t118 * t142 + (t118 * t125 - t139) * t115;
	t155 = t100 * t101 * t91;
	t99 = t103 ^ 2;
	t89 = t99 * t101 + 0.1e1;
	t90 = -t104 * t118 + t114 * t139;
	t158 = (-t90 * t152 - t99 * t155) / t89 ^ 2;
	t84 = 0.1e1 / t86;
	t156 = t84 * t93;
	t144 = -0.2e1 * t158;
	t137 = -t100 * t114 + t115 * t152;
	t87 = 0.1e1 / t89;
	t83 = 0.2e1 * (t110 - t138 * t116 / t112 ^ 2) / t128 * t162 * t163;
	t78 = t144 + 0.2e1 * (-t101 * t87 * t90 + (-t101 * t158 - t87 * t155) * t103) * t103;
	t1 = [0, t83, 0, 0, 0; 0, ((-0.2e1 * t92 * t159 + (-qJD(2) * t82 - t81) * t156) * t128 + (t93 * t159 * t160 + (t83 * t93 * t143 + t157 * t160 - qJD(2) * t92 - (t153 * qJD(2) + t140 * t96) * t106 * t154) * t84 - (t106 * t83 + (qJD(2) + (-0.2e1 * t125 + t97) * t96) * t107) * t128 * t156) * t127) * t126, 0, 0, 0; 0, (t137 * t87 * t145 + (t137 * t144 + ((-t100 * t118 - 0.2e1 * t103 * t155) * t115 + (-t115 * t90 + (-t103 * t118 + t91) * t114) * t101) * t87) * t127) * t126, t78, t78, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:34:31
	% EndTime: 2019-10-24 10:34:32
	% DurationCPUTime: 0.77s
	% Computational Cost: add. (5005->83), mult. (5935->196), div. (1261->15), fcn. (7597->9), ass. (0->94)
	t166 = qJ(3) + qJ(4);
	t159 = cos(t166);
	t168 = cos(pkin(8));
	t158 = sin(t166);
	t167 = sin(pkin(8));
	t170 = cos(qJ(2));
	t205 = t167 * t170;
	t189 = t158 * t205;
	t145 = t168 * t159 + t189;
	t169 = sin(qJ(2));
	t202 = t169 * t158;
	t135 = atan2(-t145, t202);
	t131 = sin(t135);
	t132 = cos(t135);
	t140 = t145 ^ 2;
	t156 = 0.1e1 / t158 ^ 2;
	t164 = 0.1e1 / t169 ^ 2;
	t209 = t156 * t164;
	t136 = t140 * t209 + 0.1e1;
	t133 = 0.1e1 / t136;
	t155 = 0.1e1 / t158;
	t192 = t155 * t164 * t170;
	t182 = t145 * t192 + t167;
	t122 = t182 * t133;
	t201 = -t122 + t167;
	t224 = t201 * t169 * t131 + t132 * t170;
	t203 = t168 * t170;
	t149 = t167 * t158 + t159 * t203;
	t143 = 0.1e1 / t149 ^ 2;
	t160 = t168 ^ 2;
	t162 = t169 ^ 2;
	t190 = t160 * t162 * t143;
	t139 = 0.1e1 + t190;
	t199 = qJD(2) * t170;
	t183 = t169 * t199;
	t161 = qJD(3) + qJD(4);
	t200 = qJD(2) * t169;
	t184 = t168 * t200;
	t186 = t158 * t203;
	t130 = -t161 * t186 + (t161 * t167 - t184) * t159;
	t142 = 0.1e1 / t149;
	t216 = t130 * t142 * t143;
	t194 = t162 * t216;
	t223 = (t143 * t183 - t194) * t160 / t139 ^ 2;
	t163 = 0.1e1 / t169;
	t188 = t159 * t205;
	t147 = -t168 * t158 + t188;
	t181 = t145 * t156 * t159 - t147 * t155;
	t222 = t163 * t181;
	t180 = t161 * t168 + t167 * t200;
	t127 = t158 * t180 - t161 * t188;
	t157 = t155 * t156;
	t165 = t163 / t162;
	t208 = t159 * t161;
	t191 = t164 * t208;
	t193 = t145 * t209;
	t220 = -0.2e1 * (-t127 * t193 + (-t156 * t165 * t199 - t157 * t191) * t140) / t136 ^ 2;
	t215 = t131 * t145;
	t126 = t132 * t202 - t215;
	t123 = 0.1e1 / t126;
	t124 = 0.1e1 / t126 ^ 2;
	t148 = -t167 * t159 + t186;
	t219 = 0.2e1 * t148;
	t218 = t124 * t148;
	t129 = -t149 * t161 + t158 * t184;
	t217 = t129 * t124;
	t214 = t131 * t148;
	t212 = t132 * t145;
	t211 = t132 * t148;
	t207 = t159 * t169;
	t204 = t168 * t123;
	t179 = t158 * t199 + t161 * t207;
	t115 = (t127 * t163 * t155 + t179 * t193) * t133;
	t178 = -t115 * t145 + t179;
	t111 = (-t115 * t202 + t127) * t131 + t178 * t132;
	t141 = t148 ^ 2;
	t120 = t141 * t124 + 0.1e1;
	t125 = t123 * t124;
	t198 = 0.2e1 * (-t141 * t125 * t111 - t148 * t217) / t120 ^ 2;
	t197 = t125 * t219;
	t196 = t169 * t219;
	t195 = t124 * t214;
	t187 = t169 * t204;
	t137 = 0.1e1 / t139;
	t128 = t159 * t180 + t161 * t189;
	t118 = 0.1e1 / t120;
	t117 = t133 * t222;
	t114 = -t122 * t212 + t224 * t158;
	t113 = t132 * t207 - t131 * t147 + (-t131 * t202 - t212) * t117;
	t112 = (t143 * t196 * t223 + (t196 * t216 + (t129 * t169 - t148 * t199) * t143) * t137) * t168;
	t110 = t182 * t220 + (-t127 * t192 + (-t156 * t170 * t191 + (-0.2e1 * t165 * t170 ^ 2 - t163) * t155 * qJD(2)) * t145) * t133;
	t108 = t220 * t222 + (-t181 * t164 * t199 + ((-t145 * t161 + t128) * t155 + (-0.2e1 * t145 * t157 * t208 + (t147 * t161 - t127) * t156) * t159) * t163) * t133;
	t107 = (t113 * t218 - t123 * t149) * t198 + (t113 * t217 + t130 * t123 + (t113 * t197 - t149 * t124) * t111 - (t159 * t199 - t161 * t202 - t108 * t145 + t117 * t127 + (-t117 * t202 - t147) * t115) * t124 * t211 - (t128 + (-t108 * t158 - t115 * t159) * t169 - t178 * t117) * t195) * t118;
	t1 = [0, t110, t108, t108, 0; 0, t114 * t198 * t218 + (-(-t110 * t212 + (t115 * t215 + t127 * t132) * t122) * t218 + (-t224 * t218 - t187) * t208 + (t111 * t197 + t217) * t114) * t118 + (t187 * t198 + ((-qJD(2) * t204 - (t201 * qJD(2) - t115) * t195) * t170 + (t110 * t214 + t168 * t111 - (t201 * t115 - qJD(2)) * t211) * t124 * t169) * t118) * t158, t107, t107, 0; 0, 0.2e1 * (t142 * t203 + t159 * t190) * t223 + (0.2e1 * t159 * t160 * t194 + t142 * t184 + (t130 * t203 + (t158 * t161 * t162 - 0.2e1 * t159 * t183) * t160) * t143) * t137, t112, t112, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end