% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRP3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_jacobiaD_rot_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:09
% EndTime: 2019-02-26 21:26:10
% DurationCPUTime: 0.68s
% Computational Cost: add. (813->92), mult. (2519->208), div. (480->12), fcn. (2968->9), ass. (0->95)
t132 = sin(qJ(1));
t173 = qJD(2) * t132;
t126 = t132 ^ 2;
t131 = sin(qJ(2));
t124 = 0.1e1 / t131 ^ 2;
t134 = cos(qJ(2));
t128 = t134 ^ 2;
t184 = t124 * t128;
t119 = t126 * t184 + 0.1e1;
t117 = 0.1e1 / t119;
t123 = 0.1e1 / t131;
t160 = t124 * t173;
t135 = cos(qJ(1));
t175 = qJD(1) * t135;
t161 = t134 * t175;
t91 = ((t131 * t173 - t161) * t123 + t128 * t160) * t117;
t151 = -t91 + t173;
t179 = t132 * t134;
t116 = atan2(-t179, t131);
t115 = cos(t116);
t114 = sin(t116);
t165 = t114 * t179;
t100 = t115 * t131 - t165;
t97 = 0.1e1 / t100;
t133 = cos(qJ(5));
t178 = t133 * t135;
t163 = t131 * t178;
t130 = sin(qJ(5));
t181 = t132 * t130;
t113 = t163 - t181;
t107 = 0.1e1 / t113;
t98 = 0.1e1 / t100 ^ 2;
t108 = 0.1e1 / t113 ^ 2;
t197 = t117 - 0.1e1;
t152 = -t132 * t91 + qJD(2);
t186 = t115 * t134;
t86 = t152 * t186 + (t151 * t131 - t161) * t114;
t196 = t86 * t97 * t98;
t129 = t135 ^ 2;
t96 = t128 * t129 * t98 + 0.1e1;
t94 = 0.1e1 / t96;
t195 = t94 * t98;
t194 = t97 * t94;
t180 = t132 * t133;
t182 = t131 * t135;
t112 = t130 * t182 + t180;
t106 = t112 ^ 2;
t105 = t106 * t108 + 0.1e1;
t189 = t108 * t112;
t149 = qJD(1) * t131 + qJD(5);
t150 = qJD(5) * t131 + qJD(1);
t171 = qJD(2) * t135;
t159 = t134 * t171;
t183 = t130 * t135;
t93 = -t150 * t183 + (-t149 * t132 + t159) * t133;
t192 = t107 * t108 * t93;
t92 = -qJD(5) * t163 - t130 * t159 - t133 * t175 + t149 * t181;
t193 = 0.1e1 / t105 ^ 2 * (-t106 * t192 - t92 * t189);
t191 = t135 * t98;
t190 = t107 * t130;
t188 = t112 * t133;
t187 = t114 * t132;
t185 = t123 * t128;
t177 = qJD(1) * t132;
t176 = qJD(1) * t134;
t174 = qJD(2) * t131;
t172 = qJD(2) * t134;
t162 = t132 * t175;
t166 = t98 * t174;
t170 = 0.2e1 * (-t129 * t134 * t166 + (-t129 * t196 - t98 * t162) * t128) / t96 ^ 2;
t169 = 0.2e1 * t196;
t168 = 0.2e1 * t193;
t127 = t134 * t128;
t144 = qJD(2) * (-t124 * t127 - t134) * t123;
t167 = 0.2e1 * (t126 * t144 + t162 * t184) / t119 ^ 2;
t164 = t117 * t185;
t158 = t97 * t170;
t157 = t98 * t170;
t156 = 0.1e1 + t184;
t155 = 0.2e1 * t112 * t192;
t154 = t132 * t167;
t153 = t134 * t167;
t148 = t132 * t164;
t147 = t156 * t135;
t146 = t149 * t135;
t145 = t108 * t188 - t190;
t143 = t145 * t135;
t111 = -t131 * t180 - t183;
t110 = -t131 * t181 + t178;
t104 = t156 * t132 * t117;
t102 = 0.1e1 / t105;
t90 = (t197 * t134 * t114 + t115 * t148) * t135;
t89 = t131 * t187 + t186 + (-t114 * t131 - t115 * t179) * t104;
t87 = -t156 * t154 + (qJD(1) * t147 + 0.2e1 * t132 * t144) * t117;
t1 = [t123 * t135 * t153 + (t123 * t132 * t176 + qJD(2) * t147) * t117, t87, 0, 0, 0, 0; (t174 * t194 + (t158 + (qJD(1) * t90 + t86) * t195) * t134) * t132 + (t90 * t157 * t134 + (t90 * t166 + (t90 * t169 + ((t91 * t148 + t197 * t174 + t153) * t114 + (t154 * t185 + t134 * t91 + (t127 * t160 - (t91 - 0.2e1 * t173) * t134) * t117) * t115) * t191) * t134 + (-t97 + (-(-t126 + t129) * t115 * t164 + t197 * t165) * t98) * t176) * t94) * t135 (t177 * t194 + (t158 + (qJD(2) * t89 + t86) * t195) * t135) * t131 + (t89 * t135 * t157 + (-t97 * t171 - (-t115 * t132 * t87 + t151 * t114 + (-qJD(2) * t114 - t115 * t175 + t187 * t91) * t104) * t134 * t191 + (t135 * t169 + t98 * t177) * t89) * t94 - ((-t87 + t175) * t114 + (t151 * t104 - t152) * t115) * t182 * t195) * t134, 0, 0, 0, 0; (-t107 * t110 + t111 * t189) * t168 + (t111 * t155 - t150 * t107 * t180 + (-t132 * t172 - t146) * t190 + (-t110 * t93 + t111 * t92 + t146 * t188 - (t150 * t130 - t133 * t172) * t112 * t132) * t108) * t102, t134 * t143 * t168 + (t143 * t174 + (t145 * t177 + ((qJD(5) * t107 + t155) * t133 + (t133 * t92 + (qJD(5) * t112 - t93) * t130) * t108) * t135) * t134) * t102, 0, 0, -0.2e1 * t193 + 0.2e1 * (-t102 * t108 * t92 + (-t102 * t192 - t108 * t193) * t112) * t112, 0;];
JaD_rot  = t1;
