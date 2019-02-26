% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRP4_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_jacobiaD_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:47
% EndTime: 2019-02-26 21:26:48
% DurationCPUTime: 0.54s
% Computational Cost: add. (926->76), mult. (3228->188), div. (613->15), fcn. (4191->9), ass. (0->85)
t120 = cos(pkin(9));
t172 = 0.2e1 * t120;
t121 = sin(qJ(2));
t115 = 0.1e1 / t121;
t119 = sin(pkin(9));
t124 = cos(qJ(1));
t149 = t124 * t120;
t122 = sin(qJ(1));
t123 = cos(qJ(2));
t152 = t122 * t123;
t104 = t119 * t152 + t149;
t112 = 0.1e1 / t119;
t158 = t104 * t112;
t140 = t115 * t158;
t154 = t121 * t119;
t94 = atan2(-t104, t154);
t90 = sin(t94);
t91 = cos(t94);
t113 = 0.1e1 / t119 ^ 2;
t116 = 0.1e1 / t121 ^ 2;
t99 = t104 ^ 2;
t97 = t99 * t116 * t113 + 0.1e1;
t92 = 0.1e1 / t97;
t171 = (t91 * t140 + t90) * t92 - t90;
t86 = -t90 * t104 + t91 * t154;
t83 = 0.1e1 / t86;
t150 = t124 * t119;
t138 = t123 * t150;
t107 = -t122 * t120 + t138;
t75 = t171 * t107;
t170 = 0.2e1 * t75;
t153 = t122 * t119;
t108 = t123 * t149 + t153;
t101 = 0.1e1 / t108;
t102 = 0.1e1 / t108 ^ 2;
t84 = 0.1e1 / t86 ^ 2;
t100 = t107 ^ 2;
t145 = qJD(2) * t124;
t136 = t121 * t145;
t87 = t104 * qJD(1) + t119 * t136;
t166 = t84 * t87;
t146 = qJD(2) * t123;
t161 = t121 * t90;
t164 = t104 * t91;
t137 = t116 * t146;
t148 = qJD(1) * t122;
t89 = -qJD(2) * t121 * t153 + qJD(1) * t138 - t120 * t148;
t76 = (t104 * t137 - t115 * t89) * t92 * t112;
t73 = -t76 * t164 - t90 * t89 + (t91 * t146 - t76 * t161) * t119;
t168 = t73 * t83 * t84;
t80 = t100 * t84 + 0.1e1;
t169 = (-t100 * t168 - t107 * t166) / t80 ^ 2;
t114 = t121 ^ 2;
t117 = t115 / t114;
t167 = (t104 * t116 * t89 - t117 * t99 * t146) * t113 / t97 ^ 2;
t106 = -t120 * t152 + t150;
t133 = t120 * t136;
t88 = t106 * qJD(1) - t133;
t165 = t101 * t102 * t88;
t163 = t107 * t91;
t162 = t115 * t92;
t160 = t90 * t107;
t155 = t116 * t123;
t82 = (t155 * t158 + t122) * t92;
t159 = t122 - t82;
t157 = t106 * t124;
t118 = t124 ^ 2;
t156 = t114 * t118;
t151 = t123 * t124;
t147 = qJD(1) * t124;
t130 = t114 * t122 * t147 - t118 * t121 * t146;
t134 = t156 * t165;
t139 = t102 * t156;
t98 = 0.1e1 + t139;
t144 = 0.2e1 * (-t130 * t102 - t134) / t98 ^ 2;
t143 = -0.2e1 * t167;
t142 = t84 * t169;
t141 = t84 * t160;
t135 = 0.2e1 * t83 * t169;
t131 = 0.2e1 * t115 * t167 + t92 * t137;
t95 = 0.1e1 / t98;
t78 = 0.1e1 / t80;
t74 = -t82 * t164 + (t123 * t91 + t159 * t161) * t119;
t72 = t122 * t143 + t92 * t147 + (t89 * t92 * t155 + (t143 * t155 + (-0.2e1 * t117 * t123 ^ 2 - t115) * t92 * qJD(2)) * t104) * t112;
t1 = [(t131 * t107 + t87 * t162) * t112, t72, 0, 0, 0, 0; t104 * t135 + (-t89 * t83 + (t104 * t73 + t75 * t87) * t84) * t78 + (t142 * t170 + (t168 * t170 - (-t76 * t92 * t140 + t143) * t141 - ((t92 - 0.1e1) * t76 + (-t131 * t104 + t89 * t162) * t112) * t84 * t163 + t171 * t166) * t78) * t107, t74 * t78 * t166 + (-(-t82 * t91 * t89 + (t76 * t82 * t90 - t72 * t91) * t104) * t84 * t78 + 0.2e1 * (t78 * t168 + t142) * t74) * t107 + (t124 * t135 * t121 + ((-t83 * t145 - (t159 * qJD(2) - t76) * t141) * t123 + (t83 * t148 + (t124 * t73 - (-t72 + t147) * t160 - (t159 * t76 - qJD(2)) * t163) * t84) * t121) * t78) * t119, 0, 0, 0, 0; (-t101 * t122 - t102 * t157) * t121 * t144 + (-0.2e1 * t121 * t157 * t165 + (t121 * t147 + t122 * t146) * t101 + (((-t88 + t133) * t122 - t108 * t147) * t121 + (-t121 * t148 + t123 * t145) * t106) * t102) * t95 (t101 * t151 + t120 * t139) * t144 + (t134 * t172 + (t123 * t148 + t136) * t101 + (t130 * t172 + t88 * t151) * t102) * t95, 0, 0, 0, 0;];
JaD_rot  = t1;
