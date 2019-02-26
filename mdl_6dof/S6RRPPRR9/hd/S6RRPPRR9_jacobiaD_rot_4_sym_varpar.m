% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR9
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
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR9_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:12
% EndTime: 2019-02-26 21:33:13
% DurationCPUTime: 0.51s
% Computational Cost: add. (969->70), mult. (3196->173), div. (656->14), fcn. (4222->9), ass. (0->73)
t120 = sin(qJ(2));
t121 = sin(qJ(1));
t122 = cos(qJ(2));
t123 = cos(qJ(1));
t151 = cos(pkin(6));
t136 = t123 * t151;
t103 = t120 * t136 + t121 * t122;
t119 = sin(pkin(6));
t111 = 0.1e1 / t119;
t113 = 0.1e1 / t120;
t139 = t103 * t111 * t113;
t148 = t119 * t120;
t96 = atan2(-t103, t148);
t90 = sin(t96);
t91 = cos(t96);
t112 = 0.1e1 / t119 ^ 2;
t114 = 0.1e1 / t120 ^ 2;
t99 = t103 ^ 2;
t97 = t112 * t114 * t99 + 0.1e1;
t92 = 0.1e1 / t97;
t162 = (t91 * t139 + t90) * t92 - t90;
t85 = -t103 * t90 + t91 * t148;
t82 = 0.1e1 / t85;
t116 = 0.1e1 / t121;
t117 = 0.1e1 / t121 ^ 2;
t83 = 0.1e1 / t85 ^ 2;
t144 = qJD(2) * t122;
t153 = t120 * t90;
t158 = t103 * t91;
t138 = t114 * t144;
t155 = t111 * t92;
t137 = t121 * t151;
t135 = t120 * t137;
t146 = t123 * t122;
t147 = t121 * t120;
t89 = -qJD(1) * t135 - qJD(2) * t147 + (qJD(2) * t151 + qJD(1)) * t146;
t75 = (t103 * t138 - t113 * t89) * t155;
t72 = -t75 * t158 - t90 * t89 + (t91 * t144 - t75 * t153) * t119;
t161 = t72 * t82 * t83;
t102 = -t122 * t136 + t147;
t150 = t114 * t122;
t133 = t102 * t113 + t103 * t150;
t76 = t133 * t155;
t160 = t75 * t76;
t115 = t113 * t114;
t159 = (t103 * t114 * t89 - t115 * t99 * t144) * t112 / t97 ^ 2;
t132 = t135 - t146;
t157 = t132 * t83;
t156 = t132 * t91;
t154 = t113 * t92;
t152 = t90 * t132;
t149 = t117 * t123;
t145 = qJD(1) * t123;
t101 = t132 ^ 2;
t79 = t101 * t83 + 0.1e1;
t131 = t123 * t120 + t122 * t137;
t87 = t103 * qJD(1) + t131 * qJD(2);
t143 = 0.2e1 * (-t101 * t161 + t87 * t157) / t79 ^ 2;
t142 = 0.2e1 * t161;
t141 = -0.2e1 * t159;
t100 = t131 ^ 2;
t118 = t116 * t117;
t86 = t102 * qJD(1) + t132 * qJD(2);
t98 = t100 * t112 * t117 + 0.1e1;
t140 = 0.2e1 * (-t100 * t118 * t145 - t117 * t131 * t86) * t112 / t98 ^ 2;
t130 = 0.2e1 * t113 * t159 + t92 * t138;
t94 = 0.1e1 / t98;
t88 = t131 * qJD(1) + t103 * qJD(2);
t77 = 0.1e1 / t79;
t74 = t162 * t132;
t73 = -t76 * t158 + t90 * t102 + (t122 * t91 - t76 * t153) * t119;
t71 = (t133 * t141 + (t89 * t150 + t113 * t88 + (-t102 * t150 + (-0.2e1 * t115 * t122 ^ 2 - t113) * t103) * qJD(2)) * t92) * t111;
t1 = [(-t130 * t132 + t87 * t154) * t111, t71, 0, 0, 0, 0; t103 * t82 * t143 + (-t89 * t82 + (t103 * t72 - t74 * t87) * t83) * t77 - (-t74 * t142 * t77 + (-t74 * t143 + ((-t75 * t92 * t139 + t141) * t152 + ((t92 - 0.1e1) * t75 + (-t130 * t103 + t89 * t154) * t111) * t156 + t162 * t87) * t77) * t83) * t132 (t131 * t82 - t73 * t157) * t143 + (-t73 * t132 * t142 + t86 * t82 + (t131 * t72 + t73 * t87 + (t103 * t160 + t88) * t152 + (t102 * t75 - t103 * t71 - t76 * t89) * t156) * t83 + ((-qJD(2) * t76 - t75) * t90 * t122 + (-t71 * t90 + (-qJD(2) - t160) * t91) * t120) * t119 * t157) * t77, 0, 0, 0, 0; ((-t102 * t116 - t131 * t149) * t140 + (-t86 * t149 + t116 * t88 + (-t102 * t149 - (0.2e1 * t118 * t123 ^ 2 + t116) * t131) * qJD(1)) * t94) * t111 (t116 * t87 * t94 - (t117 * t94 * t145 + t116 * t140) * t132) * t111, 0, 0, 0, 0;];
JaD_rot  = t1;
