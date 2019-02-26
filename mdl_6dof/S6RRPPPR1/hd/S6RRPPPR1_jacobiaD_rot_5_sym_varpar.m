% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:59
% EndTime: 2019-02-26 21:21:59
% DurationCPUTime: 0.56s
% Computational Cost: add. (2051->79), mult. (3228->188), div. (613->15), fcn. (4191->9), ass. (0->85)
t133 = cos(pkin(10));
t181 = 0.2e1 * t133;
t130 = qJ(2) + pkin(9);
t127 = cos(t130);
t132 = sin(pkin(10));
t134 = sin(qJ(1));
t162 = t134 * t132;
t135 = cos(qJ(1));
t163 = t133 * t135;
t114 = t127 * t162 + t163;
t126 = sin(t130);
t165 = t126 * t132;
t104 = atan2(-t114, t165);
t101 = cos(t104);
t100 = sin(t104);
t171 = t100 * t114;
t96 = t101 * t165 - t171;
t93 = 0.1e1 / t96;
t109 = t114 ^ 2;
t124 = 0.1e1 / t126 ^ 2;
t129 = 0.1e1 / t132 ^ 2;
t107 = t109 * t124 * t129 + 0.1e1;
t102 = 0.1e1 / t107;
t164 = t132 * t135;
t150 = t127 * t164;
t161 = t134 * t133;
t117 = t150 - t161;
t123 = 0.1e1 / t126;
t128 = 0.1e1 / t132;
t169 = t114 * t128;
t151 = t123 * t169;
t85 = (-t100 + (t101 * t151 + t100) * t102) * t117;
t180 = 0.2e1 * t85;
t118 = t127 * t163 + t162;
t111 = 0.1e1 / t118;
t112 = 0.1e1 / t118 ^ 2;
t94 = 0.1e1 / t96 ^ 2;
t110 = t117 ^ 2;
t156 = qJD(2) * t135;
t147 = t126 * t156;
t97 = t114 * qJD(1) + t132 * t147;
t176 = t97 * t94;
t158 = qJD(2) * t127;
t149 = t124 * t158;
t157 = qJD(2) * t134;
t160 = qJD(1) * t134;
t99 = qJD(1) * t150 - t133 * t160 - t157 * t165;
t141 = t114 * t149 - t123 * t99;
t86 = t141 * t128 * t102;
t83 = (-t114 * t86 + t132 * t158) * t101 + (-t86 * t165 - t99) * t100;
t178 = t83 * t93 * t94;
t90 = t110 * t94 + 0.1e1;
t179 = (-t110 * t178 - t117 * t176) / t90 ^ 2;
t88 = 0.1e1 / t90;
t177 = t88 * t94;
t122 = t126 ^ 2;
t125 = t123 / t122;
t175 = 0.1e1 / t107 ^ 2 * (-t109 * t125 * t158 + t114 * t124 * t99) * t129;
t116 = -t127 * t161 + t164;
t143 = t133 * t147;
t98 = t116 * qJD(1) - t143;
t174 = t111 * t112 * t98;
t173 = t123 * t97;
t166 = t124 * t127;
t142 = t166 * t169 + t134;
t92 = t142 * t102;
t172 = t134 - t92;
t170 = t100 * t117;
t168 = t116 * t135;
t131 = t135 ^ 2;
t167 = t122 * t131;
t159 = qJD(1) * t135;
t155 = -0.2e1 * t175;
t152 = t112 * t167;
t108 = 0.1e1 + t152;
t144 = t122 * t134 * t159;
t145 = t167 * t174;
t148 = qJD(2) * t126 * t131;
t154 = 0.2e1 / t108 ^ 2 * (-t145 + (t127 * t148 - t144) * t112);
t153 = t94 * t179;
t146 = 0.2e1 * t93 * t179;
t105 = 0.1e1 / t108;
t84 = -t101 * t114 * t92 + (t172 * t126 * t100 + t101 * t127) * t132;
t82 = t142 * t155 + (t159 + (t99 * t166 + (-0.2e1 * t125 * t127 ^ 2 - t123) * t114 * qJD(2)) * t128) * t102;
t1 = [(0.2e1 * t117 * t123 * t175 + (t117 * t149 + t173) * t102) * t128, t82, 0, 0, 0, 0; t114 * t146 + (-t99 * t93 + (t114 * t83 + t85 * t97) * t94) * t88 + (t153 * t180 + (t178 * t180 + ((t102 * t97 - t97 - (-t102 * t86 * t151 + t155) * t117) * t100 + (-(t151 * t155 - t86) * t117 + (-t117 * t86 + (t114 * t173 + t141 * t117) * t128) * t102) * t101) * t94) * t88) * t117, t84 * t88 * t176 + (-(t92 * t86 * t171 + (-t114 * t82 - t92 * t99) * t101) * t177 + 0.2e1 * (t88 * t178 + t153) * t84) * t117 + ((-t93 * t88 * t156 - (t172 * qJD(2) - t86) * t170 * t177) * t127 + (t135 * t146 + (t93 * t160 + (t135 * t83 - (-t82 + t159) * t170 - (t172 * t86 - qJD(2)) * t117 * t101) * t94) * t88) * t126) * t132, 0, 0, 0, 0; (-t111 * t134 - t112 * t168) * t126 * t154 + (-0.2e1 * t126 * t168 * t174 + (t126 * t159 + t127 * t157) * t111 + (((-t98 + t143) * t134 - t118 * t159) * t126 + (-t126 * t160 + t127 * t156) * t116) * t112) * t105 (t111 * t127 * t135 + t133 * t152) * t154 + (t145 * t181 + (t127 * t160 + t147) * t111 + (t144 * t181 + (-0.2e1 * t133 * t148 + t135 * t98) * t127) * t112) * t105, 0, 0, 0, 0;];
JaD_rot  = t1;
