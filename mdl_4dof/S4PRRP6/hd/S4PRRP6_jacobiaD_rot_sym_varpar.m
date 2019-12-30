% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4PRRP6
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
%   Wie in S4PRRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4PRRP6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRP6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_jacobiaD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:22:13
	% EndTime: 2019-12-29 12:22:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:22:07
	% EndTime: 2019-12-29 12:22:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:22:07
	% EndTime: 2019-12-29 12:22:07
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:22:13
	% EndTime: 2019-12-29 12:22:13
	% DurationCPUTime: 0.61s
	% Computational Cost: add. (429->41), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->59)
	t97 = sin(qJ(2));
	t90 = t97 ^ 2;
	t99 = cos(qJ(2));
	t123 = t90 / t99 ^ 2;
	t111 = 0.1e1 + t123;
	t94 = sin(pkin(6));
	t135 = t111 * t94;
	t134 = qJD(2) * t97;
	t95 = cos(pkin(6));
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
	t1 = [0, t56, 0, 0; 0, ((-0.2e1 * t60 * t131 + (-qJD(2) * t54 - t53) * t128) * t99 + (t61 * t131 * t132 + (t56 * t61 * t115 + t130 * t132 - qJD(2) * t60 - (-t119 * qJD(2) + t112 * t64) * t78 * t127) * t57 - (t56 * t78 + (-t112 * qJD(2) + t119 * t64) * t79) * t99 * t128) * t97) * t95, 0, 0; 0, (t110 * t97 * t116 + (t110 * t118 + ((-qJD(3) * t73 - 0.2e1 * t76 * t126) * t98 + (-t69 * t98 + (t70 - t117) * t96) * t74) * t97) * t65) * t95, t116 + 0.2e1 * (-t65 * t69 * t74 + (-t65 * t126 - t74 * t129) * t76) * t76, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:22:13
	% EndTime: 2019-12-29 12:22:14
	% DurationCPUTime: 1.42s
	% Computational Cost: add. (1126->84), mult. (3936->208), div. (821->15), fcn. (5024->9), ass. (0->85)
	t129 = sin(qJ(2));
	t131 = cos(qJ(2));
	t126 = sin(pkin(6));
	t127 = cos(pkin(6));
	t130 = cos(qJ(3));
	t128 = sin(qJ(3));
	t159 = t128 * t131;
	t108 = t126 * t159 + t127 * t130;
	t119 = 0.1e1 / t128;
	t124 = 0.1e1 / t129 ^ 2;
	t161 = t124 * t131;
	t146 = t119 * t161;
	t103 = t108 ^ 2;
	t120 = 0.1e1 / t128 ^ 2;
	t163 = t120 * t124;
	t102 = t103 * t163 + 0.1e1;
	t98 = 0.1e1 / t102;
	t85 = (t108 * t146 + t126) * t98;
	t165 = t126 - t85;
	t158 = t129 * t128;
	t100 = atan2(-t108, t158);
	t94 = sin(t100);
	t95 = cos(t100);
	t177 = t165 * t94 * t129 + t131 * t95;
	t157 = t130 * t131;
	t110 = t126 * t157 - t127 * t128;
	t162 = t120 * t130;
	t140 = t108 * t162 - t110 * t119;
	t176 = t140 * t98;
	t171 = t108 * t94;
	t89 = t95 * t158 - t171;
	t86 = 0.1e1 / t89;
	t112 = t126 * t128 + t127 * t157;
	t105 = 0.1e1 / t112;
	t123 = 0.1e1 / t129;
	t106 = 0.1e1 / t112 ^ 2;
	t87 = 0.1e1 / t89 ^ 2;
	t155 = qJD(2) * t131;
	t144 = t128 * t155;
	t153 = qJD(3) * t130;
	t139 = t129 * t153 + t144;
	t148 = t108 * t163;
	t141 = t131 * t153;
	t156 = qJD(2) * t129;
	t90 = -t126 * t141 + (qJD(3) * t127 + t126 * t156) * t128;
	t78 = (t90 * t123 * t119 + t139 * t148) * t98;
	t172 = t108 * t78;
	t75 = (-t78 * t158 + t90) * t94 + (t139 - t172) * t95;
	t88 = t86 * t87;
	t175 = t75 * t88;
	t145 = t127 * t156;
	t154 = qJD(3) * t128;
	t92 = -t126 * t154 - t127 * t141 + t128 * t145;
	t174 = t92 * t87;
	t111 = -t126 * t130 + t127 * t159;
	t143 = t130 * t156;
	t93 = -t111 * qJD(3) - t127 * t143;
	t173 = t105 * t106 * t93;
	t170 = t108 * t95;
	t169 = t111 * t87;
	t118 = t127 ^ 2;
	t122 = t129 ^ 2;
	t164 = t118 * t122;
	t101 = t106 * t164 + 0.1e1;
	t96 = 0.1e1 / t101;
	t167 = t129 * t96;
	t160 = t127 * t131;
	t104 = t111 ^ 2;
	t83 = t104 * t87 + 0.1e1;
	t152 = 0.2e1 * (-t104 * t175 - t92 * t169) / t83 ^ 2;
	t121 = t119 * t120;
	t125 = t123 / t122;
	t151 = -0.2e1 * (-t90 * t148 + (-t120 * t125 * t155 - t121 * t124 * t153) * t103) / t102 ^ 2;
	t150 = 0.2e1 * (t106 * t129 * t155 - t122 * t173) * t118 / t101 ^ 2;
	t149 = 0.2e1 * t173;
	t147 = t130 * t164;
	t142 = t130 * t155;
	t91 = qJD(3) * t108 + t126 * t143;
	t81 = 0.1e1 / t83;
	t80 = t123 * t176;
	t77 = t177 * t128 - t85 * t170;
	t76 = -t80 * t170 - t110 * t94 + (-t128 * t80 * t94 + t130 * t95) * t129;
	t74 = -t90 * t98 * t146 + t126 * t151 + (-t98 * t141 * t163 + (t151 * t161 + (-0.2e1 * t125 * t131 ^ 2 - t123) * t98 * qJD(2)) * t119) * t108;
	t72 = -t124 * t155 * t176 + (t140 * t151 + (-t90 * t162 + t119 * t91 + (t110 * t162 + (-0.2e1 * t121 * t130 ^ 2 - t119) * t108) * qJD(3)) * t98) * t123;
	t1 = [0, t74, t72, 0; 0, t77 * t152 * t169 + (t77 * t174 + (-(-t74 * t170 + (t171 * t78 + t90 * t95) * t85) * t87 + 0.2e1 * t77 * t175) * t111 + (-t127 * t129 * t86 - t177 * t169) * t153) * t81 + (-((t165 * qJD(2) - t78) * t94 * t131 + (-t74 * t94 + (t165 * t78 - qJD(2)) * t95) * t129) * t81 * t169 + (-t86 * t81 * t155 + (t75 * t81 * t87 + t86 * t152) * t129) * t127) * t128, (-t112 * t86 + t76 * t169) * t152 + (t76 * t174 + t93 * t86 + (0.2e1 * t76 * t111 * t88 - t112 * t87) * t75 + (-(t142 - t108 * t72 - t110 * t78 + t80 * t90 + (-t78 * t80 - qJD(3)) * t158) * t95 - (t91 + (-t144 + t172) * t80 + (-t128 * t72 + (-qJD(3) * t80 - t78) * t130) * t129) * t94) * t169) * t81, 0; 0, (t105 * t160 + t106 * t147) * t150 + (t147 * t149 + t105 * t145 + (t93 * t160 + (t122 * t154 - 0.2e1 * t129 * t142) * t118) * t106) * t96, (t111 * t149 * t167 + (t92 * t167 + (t129 * t150 - t96 * t155) * t111) * t106) * t127, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end