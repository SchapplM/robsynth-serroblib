% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRP2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:30:41
% EndTime: 2019-02-26 20:30:42
% DurationCPUTime: 0.74s
% Computational Cost: add. (3244->94), mult. (2519->205), div. (480->12), fcn. (2968->9), ass. (0->93)
t144 = qJ(1) + pkin(9);
t140 = sin(t144);
t134 = t140 ^ 2;
t143 = pkin(10) + qJ(4);
t139 = sin(t143);
t133 = t139 ^ 2;
t141 = cos(t143);
t136 = 0.1e1 / t141 ^ 2;
t193 = t133 * t136;
t128 = t134 * t193 + 0.1e1;
t132 = t139 * t133;
t135 = 0.1e1 / t141;
t190 = t135 * t139;
t154 = qJD(4) * (t132 * t135 * t136 + t190);
t142 = cos(t144);
t181 = qJD(1) * t142;
t191 = t133 * t140;
t159 = t181 * t191;
t197 = (t134 * t154 + t136 * t159) / t128 ^ 2;
t207 = -0.2e1 * t197;
t165 = 0.1e1 + t193;
t206 = t140 * t165;
t146 = cos(qJ(5));
t184 = t142 * t146;
t145 = sin(qJ(5));
t187 = t140 * t145;
t122 = t141 * t184 + t187;
t162 = qJD(5) * t141 - qJD(1);
t180 = qJD(4) * t139;
t205 = t145 * t162 + t146 * t180;
t188 = t140 * t139;
t125 = atan2(-t188, -t141);
t124 = cos(t125);
t123 = sin(t125);
t173 = t123 * t188;
t109 = -t124 * t141 - t173;
t106 = 0.1e1 / t109;
t116 = 0.1e1 / t122;
t107 = 0.1e1 / t109 ^ 2;
t117 = 0.1e1 / t122 ^ 2;
t126 = 0.1e1 / t128;
t204 = t126 - 0.1e1;
t172 = t126 * t133 * t135;
t160 = t140 * t172;
t99 = (t123 * t139 * t204 - t124 * t160) * t142;
t203 = t107 * t99;
t179 = qJD(4) * t140;
t169 = t136 * t179;
t170 = t139 * t181;
t178 = qJD(4) * t141;
t100 = (-(-t140 * t178 - t170) * t135 + t133 * t169) * t126;
t194 = t124 * t139;
t95 = (-t100 * t140 + qJD(4)) * t194 + (-t170 + (t100 - t179) * t141) * t123;
t202 = t106 * t107 * t95;
t155 = t141 * t187 + t184;
t168 = t145 * t180;
t104 = t155 * qJD(1) - t122 * qJD(5) + t142 * t168;
t185 = t142 * t145;
t186 = t140 * t146;
t121 = t141 * t185 - t186;
t115 = t121 ^ 2;
t114 = t115 * t117 + 0.1e1;
t196 = t117 * t121;
t161 = -qJD(1) * t141 + qJD(5);
t157 = t161 * t146;
t105 = t140 * t157 - t205 * t142;
t200 = t105 * t116 * t117;
t201 = 0.1e1 / t114 ^ 2 * (-t104 * t196 - t115 * t200);
t199 = t107 * t139;
t198 = t107 * t142;
t195 = t123 * t141;
t138 = t142 ^ 2;
t192 = t133 * t138;
t189 = t139 * t142;
t111 = t126 * t206;
t183 = t111 - t140;
t182 = qJD(1) * t140;
t103 = t107 * t192 + 0.1e1;
t177 = 0.2e1 / t103 ^ 2 * (-t192 * t202 + (t138 * t139 * t178 - t159) * t107);
t176 = 0.2e1 * t202;
t175 = -0.2e1 * t201;
t174 = t121 * t200;
t166 = t111 * t140 - 0.1e1;
t164 = t139 * t177;
t163 = t135 * t207;
t158 = t165 * t142;
t156 = -t116 * t145 + t146 * t196;
t120 = -t141 * t186 + t185;
t112 = 0.1e1 / t114;
t101 = 0.1e1 / t103;
t97 = -t140 * t195 + t194 + (-t124 * t188 + t195) * t111;
t96 = t206 * t207 + (qJD(1) * t158 + 0.2e1 * t140 * t154) * t126;
t1 = [t163 * t189 + (qJD(4) * t158 - t182 * t190) * t126, 0, 0, t96, 0, 0; (t106 * t164 + (-t106 * t178 + (qJD(1) * t99 + t95) * t199) * t101) * t140 + (t164 * t203 + (-t178 * t203 + (t99 * t176 + ((-t100 * t160 + 0.2e1 * t139 * t197 - t178 * t204) * t123 + (t163 * t191 + t100 * t139 + (t132 * t169 - (t100 - 0.2e1 * t179) * t139) * t126) * t124) * t198) * t139 + (-t106 + ((-t134 + t138) * t124 * t172 + t204 * t173) * t107) * t139 * qJD(1)) * t101) * t142, 0, 0 (-t106 * t141 + t199 * t97) * t142 * t177 + ((-t106 * t182 + (-qJD(4) * t97 - t95) * t198) * t141 + ((-qJD(4) * t106 + t176 * t97) * t142 + (t97 * t182 + (-(-t111 * t181 - t140 * t96) * t124 - (-qJD(4) * t183 + t100 * t166) * t123) * t189) * t107 - ((t96 - t181) * t123 + (-qJD(4) * t166 + t100 * t183) * t124) * t141 * t198) * t139) * t101, 0, 0; 0.2e1 * (t116 * t155 + t120 * t196) * t201 + (0.2e1 * t120 * t174 + (t120 * t104 + t155 * t105 + (-t205 * t140 - t142 * t157) * t121) * t117 + (t161 * t185 + (-t146 * t162 + t168) * t140) * t116) * t112, 0, 0, t156 * t175 * t189 + (t156 * t142 * t178 + (-t156 * t182 + ((-qJD(5) * t116 - 0.2e1 * t174) * t146 + (-t104 * t146 + (-qJD(5) * t121 + t105) * t145) * t117) * t142) * t139) * t112, t175 + 0.2e1 * (-t104 * t112 * t117 + (-t112 * t200 - t117 * t201) * t121) * t121, 0;];
JaD_rot  = t1;
