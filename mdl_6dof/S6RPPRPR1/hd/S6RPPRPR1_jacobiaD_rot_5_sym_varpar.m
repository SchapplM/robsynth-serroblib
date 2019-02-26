% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:25:29
% EndTime: 2019-02-26 20:25:30
% DurationCPUTime: 0.69s
% Computational Cost: add. (2921->87), mult. (2191->194), div. (456->12), fcn. (2616->9), ass. (0->91)
t116 = qJ(1) + pkin(9);
t112 = sin(t116);
t106 = t112 ^ 2;
t115 = pkin(10) + qJ(4);
t111 = sin(t115);
t105 = t111 ^ 2;
t113 = cos(t115);
t108 = 0.1e1 / t113 ^ 2;
t156 = t105 * t108;
t103 = t106 * t156 + 0.1e1;
t104 = t111 * t105;
t107 = 0.1e1 / t113;
t155 = t107 * t111;
t127 = qJD(4) * (t104 * t107 * t108 + t155);
t114 = cos(t116);
t147 = qJD(1) * t114;
t135 = t112 * t147;
t164 = 0.1e1 / t103 ^ 2 * (t106 * t127 + t135 * t156);
t177 = -0.2e1 * t164;
t101 = 0.1e1 / t103;
t130 = 0.1e1 + t156;
t173 = t112 * t130;
t84 = t101 * t173;
t176 = t112 * t84 - 0.1e1;
t118 = cos(pkin(11));
t146 = qJD(4) * t111;
t132 = t118 * t146;
t117 = sin(pkin(11));
t150 = t114 * t117;
t151 = t112 * t118;
t95 = -t113 * t151 + t150;
t89 = t95 * qJD(1) - t114 * t132;
t149 = t114 * t118;
t152 = t112 * t117;
t97 = t113 * t149 + t152;
t92 = 0.1e1 / t97 ^ 2;
t175 = t89 * t92;
t96 = t113 * t150 - t151;
t165 = t92 * t96;
t90 = t96 ^ 2;
t87 = t90 * t92 + 0.1e1;
t85 = 0.1e1 / t87;
t91 = 0.1e1 / t97;
t174 = (-t117 * t91 + t118 * t165) * t85;
t153 = t112 * t111;
t100 = atan2(-t153, -t113);
t98 = sin(t100);
t138 = t98 * t153;
t99 = cos(t100);
t82 = -t113 * t99 - t138;
t79 = 0.1e1 / t82;
t80 = 0.1e1 / t82 ^ 2;
t110 = t114 ^ 2;
t144 = qJD(4) * t113;
t137 = t80 * t144;
t145 = qJD(4) * t112;
t160 = t113 * t98;
t134 = t108 * t145;
t75 = (-(-t111 * t147 - t112 * t144) * t107 + t105 * t134) * t101;
t70 = (t75 - t145) * t160 + (-t98 * t147 + (-t112 * t75 + qJD(4)) * t99) * t111;
t171 = t70 * t79 * t80;
t78 = t105 * t110 * t80 + 0.1e1;
t172 = (t110 * t111 * t137 + (-t110 * t171 - t80 * t135) * t105) / t78 ^ 2;
t166 = t91 * t175;
t133 = t117 * t146;
t94 = -t113 * t152 - t149;
t88 = t94 * qJD(1) - t114 * t133;
t170 = (t88 * t165 - t90 * t166) / t87 ^ 2;
t76 = 0.1e1 / t78;
t168 = t76 * t80;
t167 = t79 * t76;
t163 = t111 * t98;
t162 = t111 * t99;
t159 = t114 * t80;
t158 = qJD(4) * t84;
t157 = t105 * t107;
t154 = t111 * t114;
t148 = qJD(1) * t112;
t143 = qJD(4) * t114;
t142 = 0.2e1 * t171;
t141 = 0.2e1 * t170;
t140 = t79 * t172;
t139 = t96 * t166;
t136 = t112 * t157;
t131 = 0.2e1 * t80 * t172;
t129 = t130 * t114;
t126 = -t99 * t136 + t163;
t74 = (t126 * t101 - t163) * t114;
t72 = (-t112 + t84) * t160 - t176 * t162;
t71 = t173 * t177 + (qJD(1) * t129 + 0.2e1 * t112 * t127) * t101;
t1 = [t107 * t154 * t177 + (qJD(4) * t129 - t148 * t155) * t101, 0, 0, t71, 0, 0; (-t144 * t167 + (0.2e1 * t140 + (qJD(1) * t74 + t70) * t168) * t111) * t112 + (t74 * t131 * t111 + (-t74 * t137 + (t74 * t142 + (t98 * t144 + t75 * t162 + 0.2e1 * t126 * t164 + ((-t75 * t136 - t144) * t98 + (t104 * t134 - (t75 - 0.2e1 * t145) * t111) * t99) * t101) * t159) * t111 + (-t79 + (-t138 + (t138 - (t106 - t110) * t99 * t157) * t101) * t80) * t111 * qJD(1)) * t76) * t114, 0, 0 (-t148 * t167 + (-0.2e1 * t140 + (-qJD(4) * t72 - t70) * t168) * t114) * t113 + (t72 * t114 * t131 + (-t79 * t143 - ((-t112 * t71 - t147 * t84) * t99 + (t176 * t75 + t145 - t158) * t98) * t80 * t154 + (t114 * t142 + t80 * t148) * t72 - ((t71 - t147) * t98 + (t75 * t84 + qJD(4) + (-t75 - t158) * t112) * t99) * t113 * t159) * t76) * t111, 0, 0; (t95 * t165 - t91 * t94) * t141 + ((-t96 * qJD(1) + t112 * t133) * t91 + 0.2e1 * t95 * t139 + (-t94 * t89 - (-t97 * qJD(1) + t112 * t132) * t96 - t95 * t88) * t92) * t85, 0, 0, t113 * t143 * t174 + (-t148 * t174 + ((t91 * t141 + t85 * t175) * t117 + (-0.2e1 * t165 * t170 + (t88 * t92 - 0.2e1 * t139) * t85) * t118) * t114) * t111, 0, 0;];
JaD_rot  = t1;
