% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:13
% EndTime: 2019-02-26 21:24:13
% DurationCPUTime: 0.54s
% Computational Cost: add. (926->75), mult. (3228->188), div. (613->15), fcn. (4191->9), ass. (0->83)
t122 = sin(qJ(2));
t116 = 0.1e1 / t122;
t121 = cos(pkin(9));
t120 = sin(pkin(9));
t125 = cos(qJ(1));
t150 = t125 * t120;
t123 = sin(qJ(1));
t124 = cos(qJ(2));
t152 = t123 * t124;
t106 = t121 * t152 - t150;
t113 = 0.1e1 / t121;
t157 = t106 * t113;
t140 = t116 * t157;
t154 = t122 * t121;
t95 = atan2(-t106, t154);
t91 = sin(t95);
t92 = cos(t95);
t100 = t106 ^ 2;
t114 = 0.1e1 / t121 ^ 2;
t117 = 0.1e1 / t122 ^ 2;
t98 = t100 * t114 * t117 + 0.1e1;
t93 = 0.1e1 / t98;
t171 = (t92 * t140 + t91) * t93 - t91;
t87 = -t106 * t91 + t92 * t154;
t84 = 0.1e1 / t87;
t151 = t124 * t125;
t109 = t123 * t120 + t121 * t151;
t76 = t171 * t109;
t170 = 0.2e1 * t76;
t153 = t123 * t121;
t108 = t124 * t150 - t153;
t101 = 0.1e1 / t108;
t102 = 0.1e1 / t108 ^ 2;
t85 = 0.1e1 / t87 ^ 2;
t104 = t109 ^ 2;
t146 = qJD(2) * t125;
t138 = t122 * t146;
t149 = qJD(1) * t123;
t133 = t124 * t149 + t138;
t148 = qJD(1) * t125;
t89 = -t120 * t148 + t133 * t121;
t166 = t85 * t89;
t147 = qJD(2) * t124;
t160 = t122 * t91;
t164 = t106 * t92;
t139 = t117 * t147;
t90 = -qJD(2) * t122 * t153 + t109 * qJD(1);
t77 = (t106 * t139 - t116 * t90) * t93 * t113;
t74 = -t77 * t164 - t90 * t91 + (t92 * t147 - t77 * t160) * t121;
t168 = t74 * t84 * t85;
t81 = t104 * t85 + 0.1e1;
t169 = (-t104 * t168 - t109 * t166) / t81 ^ 2;
t115 = t122 ^ 2;
t118 = t116 / t115;
t167 = (-t100 * t118 * t147 + t106 * t117 * t90) * t114 / t98 ^ 2;
t105 = -t120 * t152 - t121 * t125;
t135 = t120 * t138;
t88 = t105 * qJD(1) - t135;
t165 = t101 * t102 * t88;
t163 = t109 * t91;
t162 = t109 * t92;
t161 = t116 * t93;
t155 = t117 * t124;
t83 = (t155 * t157 + t123) * t93;
t159 = t123 - t83;
t158 = t105 * t125;
t119 = t125 ^ 2;
t156 = t115 * t119;
t131 = -t115 * t123 * t148 + t119 * t122 * t147;
t136 = t156 * t165;
t141 = t102 * t156;
t99 = 0.1e1 + t141;
t145 = 0.2e1 * (t131 * t102 - t136) / t99 ^ 2;
t144 = -0.2e1 * t167;
t143 = t85 * t169;
t142 = t85 * t163;
t137 = 0.2e1 * t84 * t169;
t132 = 0.2e1 * t116 * t167 + t93 * t139;
t96 = 0.1e1 / t99;
t79 = 0.1e1 / t81;
t75 = -t83 * t164 + (t124 * t92 + t159 * t160) * t121;
t73 = t123 * t144 + t93 * t148 + (t90 * t93 * t155 + (t144 * t155 + (-0.2e1 * t118 * t124 ^ 2 - t116) * t93 * qJD(2)) * t106) * t113;
t1 = [(t132 * t109 + t89 * t161) * t113, t73, 0, 0, 0, 0; t106 * t137 + (-t90 * t84 + (t106 * t74 + t76 * t89) * t85) * t79 + (t143 * t170 + (t168 * t170 - (-t77 * t93 * t140 + t144) * t142 - ((t93 - 0.1e1) * t77 + (-t132 * t106 + t90 * t161) * t113) * t85 * t162 + t171 * t166) * t79) * t109, t75 * t79 * t166 + (-(-t83 * t92 * t90 + (t77 * t83 * t91 - t73 * t92) * t106) * t85 * t79 + 0.2e1 * (t79 * t168 + t143) * t75) * t109 + (t125 * t137 * t122 + ((-t84 * t146 - (t159 * qJD(2) - t77) * t142) * t124 + (t84 * t149 + (t125 * t74 - (-t73 + t148) * t163 - (t159 * t77 - qJD(2)) * t162) * t85) * t122) * t79) * t121, 0, 0, 0, 0; (t101 * t123 + t102 * t158) * t122 * t145 + (0.2e1 * t122 * t158 * t165 + (-t122 * t148 - t123 * t147) * t101 + (((t88 - t135) * t123 + t108 * t148) * t122 + (t122 * t149 - t124 * t146) * t105) * t102) * t96 (-t101 * t151 - t120 * t141) * t145 + (-0.2e1 * t120 * t136 - t133 * t101 + (0.2e1 * t120 * t131 - t88 * t151) * t102) * t96, 0, 0, 0, 0;];
JaD_rot  = t1;
