% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR6_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:51:47
% EndTime: 2019-02-26 20:51:47
% DurationCPUTime: 0.69s
% Computational Cost: add. (2086->84), mult. (2191->187), div. (456->12), fcn. (2616->9), ass. (0->88)
t104 = pkin(10) + qJ(3);
t103 = cos(t104);
t100 = 0.1e1 / t103 ^ 2;
t102 = sin(t104);
t98 = t102 ^ 2;
t151 = t100 * t98;
t109 = sin(qJ(1));
t127 = 0.1e1 + t151;
t105 = t109 ^ 2;
t96 = t105 * t151 + 0.1e1;
t94 = 0.1e1 / t96;
t120 = t127 * t94;
t77 = t109 * t120;
t168 = t109 * t77 - 0.1e1;
t108 = cos(pkin(11));
t110 = cos(qJ(1));
t136 = qJD(3) * t110;
t125 = t102 * t136;
t142 = t109 * t108;
t107 = sin(pkin(11));
t146 = t107 * t110;
t88 = -t103 * t142 + t146;
t82 = t88 * qJD(1) - t108 * t125;
t143 = t109 * t107;
t145 = t108 * t110;
t90 = t103 * t145 + t143;
t85 = 0.1e1 / t90 ^ 2;
t167 = t82 * t85;
t166 = t102 * t151;
t89 = t103 * t146 - t142;
t153 = t85 * t89;
t83 = t89 ^ 2;
t80 = t83 * t85 + 0.1e1;
t78 = 0.1e1 / t80;
t84 = 0.1e1 / t90;
t165 = (-t107 * t84 + t108 * t153) * t78;
t144 = t109 * t102;
t93 = atan2(-t144, -t103);
t91 = sin(t93);
t130 = t91 * t144;
t92 = cos(t93);
t75 = -t103 * t92 - t130;
t72 = 0.1e1 / t75;
t99 = 0.1e1 / t103;
t73 = 0.1e1 / t75 ^ 2;
t164 = 0.2e1 * t102;
t163 = t94 - 0.1e1;
t106 = t110 ^ 2;
t139 = qJD(1) * t110;
t121 = t109 * t98 * t139;
t138 = qJD(3) * t103;
t128 = t73 * t138;
t137 = qJD(3) * t109;
t150 = t103 * t91;
t68 = (-(-t102 * t139 - t103 * t137) * t99 + t137 * t151) * t94;
t63 = (t68 - t137) * t150 + (-t91 * t139 + (-t109 * t68 + qJD(3)) * t92) * t102;
t161 = t63 * t72 * t73;
t156 = t73 * t98;
t71 = t106 * t156 + 0.1e1;
t162 = (-t73 * t121 + (t102 * t128 - t98 * t161) * t106) / t71 ^ 2;
t154 = t84 * t167;
t87 = -t103 * t143 - t145;
t81 = t87 * qJD(1) - t107 * t125;
t160 = (t81 * t153 - t83 * t154) / t80 ^ 2;
t69 = 0.1e1 / t71;
t158 = t69 * t73;
t157 = t72 * t69;
t118 = qJD(3) * (t102 + t166) * t99;
t155 = (t100 * t121 + t105 * t118) / t96 ^ 2;
t152 = t94 * t99;
t148 = t110 * t73;
t147 = qJD(3) * t77;
t141 = qJD(1) * t102;
t140 = qJD(1) * t109;
t135 = 0.2e1 * t161;
t134 = 0.2e1 * t160;
t133 = t72 * t162;
t132 = t89 * t154;
t131 = t109 * t152;
t129 = t163 * t102;
t126 = t102 * t137;
t124 = 0.2e1 * t73 * t162;
t123 = -0.2e1 * t99 * t155;
t122 = t98 * t131;
t67 = (-t92 * t122 + t91 * t129) * t110;
t65 = (-t109 + t77) * t150 - t168 * t92 * t102;
t64 = t120 * t139 + 0.2e1 * (t118 * t94 - t127 * t155) * t109;
t1 = [-t131 * t141 + (qJD(3) * t120 + t102 * t123) * t110, 0, t64, 0, 0, 0; (-t138 * t157 + (0.2e1 * t133 + (qJD(1) * t67 + t63) * t158) * t102) * t109 + (t67 * t124 * t102 + (-t67 * t128 + (t67 * t135 + ((-t68 * t122 - t163 * t138 + t155 * t164) * t91 + (-t68 * t129 + (t98 * t123 + (t164 + t166) * t94 * qJD(3)) * t109) * t92) * t148) * t102 + (-t72 - (t105 - t106) * t92 * t152 * t156 + t163 * t73 * t130) * t141) * t69) * t110, 0 (-t140 * t157 + (-0.2e1 * t133 + (-qJD(3) * t65 - t63) * t158) * t110) * t103 + (t65 * t110 * t124 + (-t72 * t136 + (t110 * t135 + t73 * t140) * t65 + (-((-t109 * t64 - t139 * t77) * t92 + (t168 * t68 + t137 - t147) * t91) * t102 - ((t64 - t139) * t91 + (t68 * t77 + qJD(3) + (-t68 - t147) * t109) * t92) * t103) * t148) * t69) * t102, 0, 0, 0; (t88 * t153 - t84 * t87) * t134 + ((-t89 * qJD(1) + t107 * t126) * t84 + 0.2e1 * t88 * t132 + (-t87 * t82 - (-t90 * qJD(1) + t108 * t126) * t89 - t88 * t81) * t85) * t78, 0, t103 * t136 * t165 + (-t140 * t165 + ((t84 * t134 + t78 * t167) * t107 + (-0.2e1 * t153 * t160 + (t81 * t85 - 0.2e1 * t132) * t78) * t108) * t110) * t102, 0, 0, 0;];
JaD_rot  = t1;
