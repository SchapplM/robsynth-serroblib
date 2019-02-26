% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function JaD_rot = S6RPRPRR6_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:51:42
% EndTime: 2019-02-26 20:51:43
% DurationCPUTime: 0.72s
% Computational Cost: add. (2618->94), mult. (2519->205), div. (480->12), fcn. (2968->9), ass. (0->94)
t145 = sin(qJ(1));
t143 = t145 ^ 2;
t142 = pkin(10) + qJ(3);
t138 = sin(t142);
t133 = t138 ^ 2;
t140 = cos(t142);
t135 = 0.1e1 / t140 ^ 2;
t192 = t133 * t135;
t128 = t143 * t192 + 0.1e1;
t132 = t138 * t133;
t134 = 0.1e1 / t140;
t189 = t134 * t138;
t155 = qJD(3) * (t132 * t134 * t135 + t189);
t146 = cos(qJ(1));
t180 = qJD(1) * t146;
t190 = t133 * t145;
t159 = t180 * t190;
t198 = (t135 * t159 + t143 * t155) / t128 ^ 2;
t207 = -0.2e1 * t198;
t165 = 0.1e1 + t192;
t206 = t145 * t165;
t141 = pkin(11) + qJ(5);
t139 = cos(t141);
t137 = sin(t141);
t185 = t145 * t137;
t186 = t140 * t146;
t122 = t139 * t186 + t185;
t184 = t145 * t138;
t125 = atan2(-t184, -t140);
t124 = cos(t125);
t123 = sin(t125);
t172 = t123 * t184;
t109 = -t124 * t140 - t172;
t106 = 0.1e1 / t109;
t116 = 0.1e1 / t122;
t107 = 0.1e1 / t109 ^ 2;
t117 = 0.1e1 / t122 ^ 2;
t126 = 0.1e1 / t128;
t205 = t126 - 0.1e1;
t171 = t126 * t133 * t134;
t160 = t145 * t171;
t99 = (t123 * t138 * t205 - t124 * t160) * t146;
t204 = t107 * t99;
t178 = qJD(3) * t145;
t168 = t135 * t178;
t169 = t138 * t180;
t100 = (-(-t140 * t178 - t169) * t134 + t133 * t168) * t126;
t193 = t124 * t138;
t95 = (-t100 * t145 + qJD(3)) * t193 + (-t169 + (t100 - t178) * t140) * t123;
t203 = t106 * t107 * t95;
t156 = t139 * t146 + t140 * t185;
t177 = qJD(3) * t146;
t167 = t138 * t177;
t101 = t156 * qJD(1) - t122 * qJD(5) + t137 * t167;
t183 = t145 * t139;
t121 = t137 * t186 - t183;
t115 = t121 ^ 2;
t114 = t115 * t117 + 0.1e1;
t196 = t117 * t121;
t161 = -qJD(1) * t140 + qJD(5);
t162 = qJD(5) * t140 - qJD(1);
t188 = t137 * t146;
t102 = -t162 * t188 + (t145 * t161 - t167) * t139;
t201 = t102 * t116 * t117;
t202 = 0.1e1 / t114 ^ 2 * (-t101 * t196 - t115 * t201);
t200 = t107 * t138;
t199 = t107 * t146;
t197 = t116 * t137;
t195 = t121 * t139;
t194 = t123 * t140;
t144 = t146 ^ 2;
t191 = t133 * t144;
t187 = t138 * t146;
t113 = t126 * t206;
t182 = t113 - t145;
t181 = qJD(1) * t145;
t179 = qJD(3) * t140;
t105 = t107 * t191 + 0.1e1;
t176 = 0.2e1 / t105 ^ 2 * (-t191 * t203 + (t138 * t144 * t179 - t159) * t107);
t175 = 0.2e1 * t203;
t174 = -0.2e1 * t202;
t173 = t121 * t201;
t166 = t113 * t145 - 0.1e1;
t164 = t138 * t176;
t163 = t134 * t207;
t158 = t165 * t146;
t157 = t117 * t195 - t197;
t154 = t138 * t178 + t146 * t161;
t120 = -t140 * t183 + t188;
t111 = 0.1e1 / t114;
t103 = 0.1e1 / t105;
t98 = -t145 * t194 + t193 + (-t124 * t184 + t194) * t113;
t96 = t206 * t207 + (qJD(1) * t158 + 0.2e1 * t145 * t155) * t126;
t1 = [t163 * t187 + (qJD(3) * t158 - t181 * t189) * t126, 0, t96, 0, 0, 0; (t106 * t164 + (-t106 * t179 + (qJD(1) * t99 + t95) * t200) * t103) * t145 + (t164 * t204 + (-t179 * t204 + (t99 * t175 + ((-t100 * t160 + 0.2e1 * t138 * t198 - t179 * t205) * t123 + (t163 * t190 + t100 * t138 + (t132 * t168 - (t100 - 0.2e1 * t178) * t138) * t126) * t124) * t199) * t138 + (-t106 + ((-t143 + t144) * t124 * t171 + t205 * t172) * t107) * t138 * qJD(1)) * t103) * t146, 0 (-t106 * t140 + t200 * t98) * t146 * t176 + ((-t106 * t181 + (-qJD(3) * t98 - t95) * t199) * t140 + ((-qJD(3) * t106 + t175 * t98) * t146 + (t98 * t181 - ((t96 - t180) * t123 + (-qJD(3) * t166 + t100 * t182) * t124) * t186 + (-(-t113 * t180 - t145 * t96) * t124 - (-qJD(3) * t182 + t100 * t166) * t123) * t187) * t107) * t138) * t103, 0, 0, 0; 0.2e1 * (t116 * t156 + t120 * t196) * t202 + (0.2e1 * t120 * t173 - t162 * t116 * t183 + t154 * t197 + (-t121 * t162 * t185 + t120 * t101 + t102 * t156 - t154 * t195) * t117) * t111, 0, t157 * t174 * t187 + (t157 * t140 * t177 + (-t157 * t181 + ((-qJD(5) * t116 - 0.2e1 * t173) * t139 + (-t101 * t139 + (-qJD(5) * t121 + t102) * t137) * t117) * t146) * t138) * t111, 0, t174 + 0.2e1 * (-t101 * t111 * t117 + (-t111 * t201 - t117 * t202) * t121) * t121, 0;];
JaD_rot  = t1;
