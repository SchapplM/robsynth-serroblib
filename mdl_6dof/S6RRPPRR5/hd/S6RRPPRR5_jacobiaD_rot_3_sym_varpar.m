% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR5_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobiaD_rot_3_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:58
% EndTime: 2019-02-26 21:30:58
% DurationCPUTime: 0.50s
% Computational Cost: add. (1133->73), mult. (3290->173), div. (633->14), fcn. (4296->9), ass. (0->78)
t123 = sin(qJ(2));
t124 = sin(qJ(1));
t125 = cos(qJ(2));
t126 = cos(qJ(1));
t157 = cos(pkin(6));
t141 = t126 * t157;
t134 = -t123 * t141 - t124 * t125;
t170 = t134 * qJD(1);
t139 = t125 * t141;
t153 = t123 * t124;
t105 = -t139 + t153;
t122 = sin(pkin(6));
t116 = 0.1e1 / t122;
t119 = 0.1e1 / t125;
t144 = t105 * t116 * t119;
t154 = t122 * t125;
t93 = atan2(-t105, -t154);
t91 = sin(t93);
t92 = cos(t93);
t100 = t105 ^ 2;
t117 = 0.1e1 / t122 ^ 2;
t120 = 0.1e1 / t125 ^ 2;
t99 = t100 * t117 * t120 + 0.1e1;
t96 = 0.1e1 / t99;
t169 = (t92 * t144 - t91) * t96 + t91;
t86 = -t105 * t91 - t92 * t154;
t83 = 0.1e1 / t86;
t142 = t124 * t157;
t140 = t123 * t142;
t152 = t126 * t125;
t109 = -t140 + t152;
t102 = 0.1e1 / t109;
t103 = 0.1e1 / t109 ^ 2;
t84 = 0.1e1 / t86 ^ 2;
t150 = qJD(2) * t123;
t159 = t125 * t91;
t165 = t105 * t92;
t143 = t120 * t150;
t162 = t116 * t96;
t135 = -t126 * t123 - t125 * t142;
t89 = -t135 * qJD(1) - t134 * qJD(2);
t76 = (t105 * t143 + t119 * t89) * t162;
t73 = -t76 * t165 - t91 * t89 + (t92 * t150 + t76 * t159) * t122;
t168 = t73 * t83 * t84;
t155 = t120 * t123;
t136 = t105 * t155 - t119 * t134;
t77 = t136 * t162;
t167 = t76 * t77;
t88 = t135 * qJD(2) + t170;
t166 = t102 * t103 * t88;
t164 = t135 * t84;
t163 = t135 * t92;
t161 = t119 * t96;
t115 = t122 ^ 2;
t118 = t124 ^ 2;
t98 = t103 * t115 * t118 + 0.1e1;
t94 = 0.1e1 / t98;
t160 = t124 * t94;
t158 = t91 * t135;
t156 = t103 * t124;
t151 = qJD(1) * t126;
t101 = t135 ^ 2;
t80 = t101 * t84 + 0.1e1;
t138 = qJD(2) * t157 + qJD(1);
t87 = -qJD(1) * t139 - qJD(2) * t152 + t138 * t153;
t149 = 0.2e1 * (-t101 * t168 + t87 * t164) / t80 ^ 2;
t148 = 0.2e1 * t168;
t147 = 0.2e1 * (-t118 * t166 + t151 * t156) * t115 / t98 ^ 2;
t121 = t119 * t120;
t146 = -0.2e1 * (t100 * t121 * t150 + t105 * t120 * t89) * t117 / t99 ^ 2;
t145 = 0.2e1 * t166;
t133 = t119 * t146 + t96 * t143;
t90 = -qJD(1) * t140 - t124 * t150 + t138 * t152;
t78 = 0.1e1 / t80;
t75 = t169 * t135;
t74 = -t77 * t165 + t91 * t134 + (t123 * t92 + t77 * t159) * t122;
t72 = (t136 * t146 + (t89 * t155 + t119 * t90 + (-t134 * t155 + (0.2e1 * t121 * t123 ^ 2 + t119) * t105) * qJD(2)) * t96) * t116;
t1 = [(-t133 * t135 - t87 * t161) * t116, t72, 0, 0, 0, 0; t105 * t83 * t149 + (-t89 * t83 + (t105 * t73 + t75 * t87) * t84) * t78 - (t75 * t148 * t78 + (t75 * t149 + ((t76 * t96 * t144 + t146) * t158 + ((t96 - 0.1e1) * t76 + (-t133 * t105 - t89 * t161) * t116) * t163 - t169 * t87) * t78) * t84) * t135 (-t109 * t83 - t74 * t164) * t149 + (-t74 * t135 * t148 + t88 * t83 + (-t109 * t73 + t74 * t87 + (t105 * t167 - t90) * t158 + (-t105 * t72 + t134 * t76 - t77 * t89) * t163) * t84 + ((-qJD(2) * t77 - t76) * t91 * t123 + (t72 * t91 + (qJD(2) + t167) * t92) * t125) * t122 * t164) * t78, 0, 0, 0, 0; ((t102 * t126 - t134 * t156) * t147 + ((qJD(1) * t102 - t134 * t145) * t124 + (-t124 * t90 + (t88 + t170) * t126) * t103) * t94) * t122 (-t135 * t145 * t160 + (t87 * t160 - (t124 * t147 - t94 * t151) * t135) * t103) * t122, 0, 0, 0, 0;];
JaD_rot  = t1;
