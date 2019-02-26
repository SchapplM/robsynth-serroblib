% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:03:28
% EndTime: 2019-02-26 22:03:29
% DurationCPUTime: 0.63s
% Computational Cost: add. (5265->87), mult. (3414->191), div. (723->12), fcn. (4027->9), ass. (0->91)
t156 = sin(qJ(1));
t211 = 0.2e1 * t156;
t150 = qJ(2) + qJ(3) + pkin(10);
t149 = cos(t150);
t155 = cos(pkin(11));
t185 = t156 * t155;
t154 = sin(pkin(11));
t157 = cos(qJ(1));
t189 = t154 * t157;
t137 = t149 * t189 - t185;
t131 = t137 ^ 2;
t186 = t156 * t154;
t188 = t155 * t157;
t138 = t149 * t188 + t186;
t133 = 0.1e1 / t138 ^ 2;
t126 = t131 * t133 + 0.1e1;
t135 = -t149 * t186 - t188;
t148 = sin(t150);
t151 = qJD(2) + qJD(3);
t190 = t151 * t157;
t175 = t148 * t190;
t127 = t135 * qJD(1) - t154 * t175;
t199 = t133 * t137;
t136 = -t149 * t185 + t189;
t128 = t136 * qJD(1) - t155 * t175;
t132 = 0.1e1 / t138;
t203 = t128 * t132 * t133;
t210 = (t127 * t199 - t131 * t203) / t126 ^ 2;
t152 = t156 ^ 2;
t144 = t148 ^ 2;
t146 = 0.1e1 / t149 ^ 2;
t197 = t144 * t146;
t142 = t152 * t197 + 0.1e1;
t140 = 0.1e1 / t142;
t145 = 0.1e1 / t149;
t183 = qJD(1) * t157;
t174 = t148 * t183;
t191 = t151 * t156;
t177 = t146 * t191;
t114 = (-(-t149 * t191 - t174) * t145 + t144 * t177) * t140;
t209 = t114 - t191;
t187 = t156 * t148;
t139 = atan2(-t187, -t149);
t130 = cos(t139);
t129 = sin(t139);
t178 = t129 * t187;
t122 = -t130 * t149 - t178;
t119 = 0.1e1 / t122;
t120 = 0.1e1 / t122 ^ 2;
t208 = t140 - 0.1e1;
t201 = t130 * t148;
t109 = (-t114 * t156 + t151) * t201 + (t209 * t149 - t174) * t129;
t207 = t109 * t119 * t120;
t143 = t148 * t144;
t194 = t145 * t148;
t165 = t151 * (t143 * t145 * t146 + t194);
t195 = t144 * t156;
t168 = t183 * t195;
t206 = (t146 * t168 + t152 * t165) / t142 ^ 2;
t205 = t120 * t148;
t204 = t120 * t157;
t202 = t129 * t156;
t200 = t132 * t154;
t198 = t144 * t145;
t153 = t157 ^ 2;
t196 = t144 * t153;
t193 = t148 * t157;
t192 = t149 * t151;
t184 = qJD(1) * t156;
t117 = t120 * t196 + 0.1e1;
t182 = 0.2e1 * (-t196 * t207 + (t148 * t153 * t192 - t168) * t120) / t117 ^ 2;
t181 = 0.2e1 * t207;
t180 = t120 * t193;
t179 = t137 * t203;
t176 = t151 * t187;
t173 = 0.1e1 + t197;
t172 = t148 * t182;
t171 = -0.2e1 * t148 * t206;
t170 = t206 * t211;
t169 = t130 * t140 * t198;
t167 = t173 * t157;
t166 = t155 * t199 - t200;
t124 = 0.1e1 / t126;
t123 = t173 * t156 * t140;
t115 = 0.1e1 / t117;
t113 = (t208 * t148 * t129 - t156 * t169) * t157;
t111 = -t149 * t202 + t201 + (t129 * t149 - t130 * t187) * t123;
t110 = -t173 * t170 + (qJD(1) * t167 + t165 * t211) * t140;
t107 = -0.2e1 * t166 * t193 * t210 + (t166 * t149 * t190 + (-0.2e1 * t179 * t188 + t184 * t200 + (t128 * t189 + (t127 * t157 - t137 * t184) * t155) * t133) * t148) * t124;
t106 = (t111 * t205 - t119 * t149) * t157 * t182 + ((-t119 * t184 + (-t111 * t151 - t109) * t204) * t149 + (-t119 * t190 - (-t110 * t130 * t156 - t209 * t129 + (t114 * t202 - t129 * t151 - t130 * t183) * t123) * t180 + (t120 * t184 + t157 * t181) * t111 - ((t110 - t183) * t129 + ((-t123 * t156 + 0.1e1) * t151 + (t123 - t156) * t114) * t130) * t149 * t204) * t148) * t115;
t1 = [t157 * t145 * t171 + (t151 * t167 - t184 * t194) * t140, t110, t110, 0, 0, 0; (t119 * t172 + (-t119 * t192 + (qJD(1) * t113 + t109) * t205) * t115) * t156 + (t120 * t172 * t113 + (-((t171 - t192 + (t114 * t145 * t195 + t192) * t140) * t129 + (t170 * t198 - t114 * t148 + (-t143 * t177 + (t114 - 0.2e1 * t191) * t148) * t140) * t130) * t180 + (-t120 * t192 + t148 * t181) * t113 + (-t119 + ((-t152 + t153) * t169 + t208 * t178) * t120) * t148 * qJD(1)) * t115) * t157, t106, t106, 0, 0, 0; 0.2e1 * (-t132 * t135 + t136 * t199) * t210 + ((-t137 * qJD(1) + t154 * t176) * t132 + 0.2e1 * t136 * t179 + (-t135 * t128 - (-t138 * qJD(1) + t155 * t176) * t137 - t136 * t127) * t133) * t124, t107, t107, 0, 0, 0;];
JaD_rot  = t1;
