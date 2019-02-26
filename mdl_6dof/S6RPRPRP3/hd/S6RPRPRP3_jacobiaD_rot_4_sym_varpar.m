% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRP3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP3_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:44:46
% EndTime: 2019-02-26 20:44:47
% DurationCPUTime: 0.61s
% Computational Cost: add. (1727->84), mult. (2191->201), div. (456->12), fcn. (2616->9), ass. (0->86)
t108 = sin(qJ(3));
t102 = t108 ^ 2;
t101 = t108 * t102;
t109 = cos(qJ(3));
t103 = 0.1e1 / t109;
t104 = 0.1e1 / t109 ^ 2;
t117 = qJD(3) * (t101 * t104 + t108) * t103;
t100 = qJ(1) + pkin(9);
t99 = cos(t100);
t142 = qJD(1) * t99;
t98 = sin(t100);
t130 = t98 * t142;
t140 = t102 * t104;
t96 = t98 ^ 2;
t95 = t96 * t140 + 0.1e1;
t151 = (t96 * t117 + t130 * t140) / t95 ^ 2;
t161 = -0.2e1 * t151;
t106 = sin(pkin(10));
t107 = cos(pkin(10));
t138 = t107 * t109;
t89 = t98 * t106 + t99 * t138;
t84 = 0.1e1 / t89 ^ 2;
t139 = t106 * t109;
t88 = -t98 * t107 + t99 * t139;
t148 = t84 * t88;
t136 = qJD(3) * t108;
t124 = t107 * t136;
t87 = t99 * t106 - t98 * t138;
t81 = t87 * qJD(1) - t99 * t124;
t159 = t81 * t84;
t83 = 0.1e1 / t89;
t149 = t83 * t159;
t82 = t88 ^ 2;
t77 = t82 * t84 + 0.1e1;
t125 = t106 * t136;
t86 = -t99 * t107 - t98 * t139;
t80 = t86 * qJD(1) - t99 * t125;
t160 = (t80 * t148 - t82 * t149) / t77 ^ 2;
t122 = 0.1e1 + t140;
t158 = t122 * t98;
t144 = t98 * t108;
t92 = atan2(-t144, -t109);
t90 = sin(t92);
t131 = t90 * t144;
t91 = cos(t92);
t76 = -t109 * t91 - t131;
t71 = 0.1e1 / t76;
t157 = -0.2e1 * t99;
t72 = 0.1e1 / t76 ^ 2;
t93 = 0.1e1 / t95;
t156 = t93 - 0.1e1;
t135 = qJD(3) * t109;
t126 = t72 * t135;
t141 = qJD(3) * t98;
t145 = t109 * t90;
t127 = t104 * t141;
t137 = qJD(1) * t108;
t67 = (-(-t98 * t135 - t99 * t137) * t103 + t102 * t127) * t93;
t62 = (t67 - t141) * t145 + (-t90 * t142 + (-t67 * t98 + qJD(3)) * t91) * t108;
t154 = t62 * t71 * t72;
t97 = t99 ^ 2;
t70 = t102 * t72 * t97 + 0.1e1;
t155 = (t108 * t97 * t126 + (-t72 * t130 - t97 * t154) * t102) / t70 ^ 2;
t153 = t67 * t90;
t152 = t72 * t99;
t79 = t93 * t158;
t150 = t79 * t98;
t147 = t79 - t98;
t146 = t108 * t99;
t143 = qJD(1) * t98;
t134 = 0.2e1 * t154;
t133 = t71 * t155;
t132 = t88 * t149;
t129 = t102 * t103 * t93;
t128 = 0.1e1 - t150;
t123 = 0.2e1 * t72 * t155;
t121 = t103 * t161;
t120 = t98 * t129;
t119 = t122 * t99;
t118 = -t106 * t83 + t107 * t148;
t74 = 0.1e1 / t77;
t68 = 0.1e1 / t70;
t66 = (t156 * t90 * t108 - t91 * t120) * t99;
t64 = t128 * t91 * t108 + t147 * t145;
t63 = t158 * t161 + (qJD(1) * t119 + 0.2e1 * t117 * t98) * t93;
t1 = [t121 * t146 + (-t103 * t98 * t137 + qJD(3) * t119) * t93, 0, t63, 0, 0, 0; (-t71 * t68 * t135 + (0.2e1 * t133 + (qJD(1) * t66 + t62) * t72 * t68) * t108) * t98 + (t66 * t123 * t108 + (-t66 * t126 + (t66 * t134 + ((0.2e1 * t108 * t151 - t67 * t120 - t156 * t135) * t90 + (t102 * t98 * t121 + t67 * t108 + (t101 * t127 - (t67 - 0.2e1 * t141) * t108) * t93) * t91) * t152) * t108 + (-t71 + (t156 * t131 - (t96 - t97) * t91 * t129) * t72) * t137) * t68) * t99, 0 (t133 * t157 + (-t71 * t143 + (-qJD(3) * t64 - t62) * t152) * t68) * t109 + (t64 * t99 * t123 + (-t99 * qJD(3) * t71 - (-t63 * t91 * t98 + t90 * t141 + t150 * t153 - t153 + (-qJD(3) * t90 - t142 * t91) * t79) * t72 * t146 + (t99 * t134 + t72 * t143) * t64 - ((t63 - t142) * t90 + (t128 * qJD(3) + t147 * t67) * t91) * t109 * t152) * t68) * t108, 0, 0, 0; 0.2e1 * (t87 * t148 - t83 * t86) * t160 + ((-t88 * qJD(1) + t98 * t125) * t83 + 0.2e1 * t87 * t132 + (-t86 * t81 - (-t89 * qJD(1) + t98 * t124) * t88 - t87 * t80) * t84) * t74, 0, t118 * t99 * t74 * t135 + (t118 * t157 * t160 + ((t83 * t143 + t99 * t159) * t106 + (t132 * t157 + (-t88 * t143 + t80 * t99) * t84) * t107) * t74) * t108, 0, 0, 0;];
JaD_rot  = t1;
