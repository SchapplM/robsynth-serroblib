% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:23:35
% EndTime: 2019-02-26 21:23:36
% DurationCPUTime: 0.51s
% Computational Cost: add. (926->75), mult. (3228->190), div. (613->15), fcn. (4191->9), ass. (0->83)
t119 = sin(pkin(9));
t121 = sin(qJ(2));
t124 = cos(qJ(1));
t120 = cos(pkin(9));
t122 = sin(qJ(1));
t151 = t122 * t120;
t106 = t119 * t124 + t121 * t151;
t123 = cos(qJ(2));
t150 = t123 * t120;
t94 = atan2(t106, t150);
t90 = sin(t94);
t91 = cos(t94);
t86 = t106 * t90 + t150 * t91;
t83 = 0.1e1 / t86;
t152 = t121 * t124;
t137 = t120 * t152;
t153 = t119 * t122;
t104 = -t137 + t153;
t115 = 0.1e1 / t123;
t112 = 0.1e1 / t120;
t157 = t106 * t112;
t138 = t115 * t157;
t103 = t106 ^ 2;
t113 = 0.1e1 / t120 ^ 2;
t116 = 0.1e1 / t123 ^ 2;
t97 = t103 * t113 * t116 + 0.1e1;
t92 = 0.1e1 / t97;
t131 = -t90 + (-t138 * t91 + t90) * t92;
t75 = t131 * t104;
t169 = 0.2e1 * t75;
t105 = t119 * t152 + t151;
t100 = 0.1e1 / t105;
t101 = 0.1e1 / t105 ^ 2;
t84 = 0.1e1 / t86 ^ 2;
t145 = qJD(2) * t124;
t135 = t123 * t145;
t88 = qJD(1) * t106 - t120 * t135;
t165 = t84 * t88;
t147 = qJD(2) * t121;
t159 = t123 * t90;
t161 = t106 * t91;
t136 = t116 * t147;
t146 = qJD(2) * t123;
t149 = qJD(1) * t122;
t87 = -qJD(1) * t137 + t119 * t149 - t146 * t151;
t76 = (t106 * t136 - t115 * t87) * t92 * t112;
t73 = t76 * t161 - t87 * t90 + (-t147 * t91 - t159 * t76) * t120;
t167 = t73 * t83 * t84;
t99 = t104 ^ 2;
t80 = t84 * t99 + 0.1e1;
t168 = (t104 * t165 - t99 * t167) / t80 ^ 2;
t114 = t123 ^ 2;
t117 = t115 / t114;
t166 = (t103 * t117 * t147 - t106 * t116 * t87) * t113 / t97 ^ 2;
t107 = t120 * t124 - t121 * t153;
t133 = t119 * t135;
t89 = qJD(1) * t107 + t133;
t164 = t100 * t101 * t89;
t163 = t104 * t90;
t162 = t104 * t91;
t160 = t115 * t92;
t154 = t116 * t121;
t82 = (t154 * t157 + t122) * t92;
t158 = t122 - t82;
t156 = t107 * t124;
t118 = t124 ^ 2;
t155 = t114 * t118;
t148 = qJD(1) * t124;
t130 = -t114 * t122 * t148 - t118 * t121 * t146;
t134 = t155 * t164;
t139 = t101 * t155;
t98 = 0.1e1 + t139;
t144 = 0.2e1 * (t101 * t130 - t134) / t98 ^ 2;
t143 = -0.2e1 * t166;
t142 = t83 * t168;
t141 = t84 * t168;
t140 = t84 * t163;
t132 = 0.2e1 * t115 * t166 - t136 * t92;
t95 = 0.1e1 / t98;
t78 = 0.1e1 / t80;
t74 = t82 * t161 + (-t121 * t91 + t158 * t159) * t120;
t72 = t122 * t143 + t92 * t148 + (-t87 * t92 * t154 + (t143 * t154 + (0.2e1 * t117 * t121 ^ 2 + t115) * t92 * qJD(2)) * t106) * t112;
t1 = [(t104 * t132 - t160 * t88) * t112, t72, 0, 0, 0, 0; -0.2e1 * t106 * t142 + (-t87 * t83 + (-t106 * t73 - t75 * t88) * t84) * t78 + (t141 * t169 + (t167 * t169 - (t138 * t76 * t92 + t143) * t140 - ((t92 - 0.1e1) * t76 + (t106 * t132 + t160 * t87) * t112) * t84 * t162 - t131 * t165) * t78) * t104, -t74 * t78 * t165 + (-(-t82 * t91 * t87 + (-t76 * t82 * t90 + t72 * t91) * t106) * t84 * t78 + 0.2e1 * (t167 * t78 + t141) * t74) * t104 + (0.2e1 * t124 * t142 * t123 + ((t83 * t145 - (-qJD(2) * t158 + t76) * t140) * t121 + (t83 * t149 + (t124 * t73 - (-t72 + t148) * t163 - (t158 * t76 - qJD(2)) * t162) * t84) * t123) * t78) * t120, 0, 0, 0, 0; (-t100 * t122 - t101 * t156) * t123 * t144 + (-0.2e1 * t123 * t156 * t164 + (-t122 * t147 + t123 * t148) * t100 + (((-t89 - t133) * t122 - t105 * t148) * t123 + (-t121 * t145 - t123 * t149) * t107) * t101) * t95 (-t100 * t152 - t119 * t139) * t144 + (-0.2e1 * t119 * t134 + (-t121 * t149 + t135) * t100 + (0.2e1 * t119 * t130 - t152 * t89) * t101) * t95, 0, 0, 0, 0;];
JaD_rot  = t1;
