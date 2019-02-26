% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:28:54
% EndTime: 2019-02-26 21:28:55
% DurationCPUTime: 0.71s
% Computational Cost: add. (2618->95), mult. (2519->203), div. (480->12), fcn. (2968->9), ass. (0->93)
t145 = qJ(2) + pkin(10);
t141 = sin(t145);
t136 = t141 ^ 2;
t143 = cos(t145);
t138 = 0.1e1 / t143 ^ 2;
t193 = t136 * t138;
t148 = sin(qJ(1));
t210 = 0.2e1 * t148;
t209 = t141 * t193;
t146 = t148 ^ 2;
t131 = t146 * t193 + 0.1e1;
t129 = 0.1e1 / t131;
t137 = 0.1e1 / t143;
t149 = cos(qJ(1));
t183 = qJD(1) * t149;
t171 = t141 * t183;
t181 = qJD(2) * t148;
t103 = (-(-t143 * t181 - t171) * t137 + t181 * t193) * t129;
t208 = t103 - t181;
t144 = pkin(11) + qJ(5);
t142 = cos(t144);
t140 = sin(t144);
t187 = t148 * t140;
t188 = t143 * t149;
t125 = t142 * t188 + t187;
t186 = t148 * t141;
t128 = atan2(-t186, -t143);
t127 = cos(t128);
t126 = sin(t128);
t174 = t126 * t186;
t112 = -t127 * t143 - t174;
t109 = 0.1e1 / t112;
t119 = 0.1e1 / t125;
t110 = 0.1e1 / t112 ^ 2;
t120 = 0.1e1 / t125 ^ 2;
t207 = -0.2e1 * t141;
t206 = t129 - 0.1e1;
t195 = t127 * t141;
t98 = (-t103 * t148 + qJD(2)) * t195 + (t208 * t143 - t171) * t126;
t205 = t109 * t110 * t98;
t159 = t142 * t149 + t143 * t187;
t180 = qJD(2) * t149;
t170 = t141 * t180;
t104 = t159 * qJD(1) - t125 * qJD(5) + t140 * t170;
t185 = t148 * t142;
t124 = t140 * t188 - t185;
t118 = t124 ^ 2;
t117 = t118 * t120 + 0.1e1;
t198 = t120 * t124;
t164 = -qJD(1) * t143 + qJD(5);
t165 = qJD(5) * t143 - qJD(1);
t190 = t140 * t149;
t105 = -t165 * t190 + (t164 * t148 - t170) * t142;
t202 = t105 * t119 * t120;
t204 = (-t104 * t198 - t118 * t202) / t117 ^ 2;
t203 = t103 * t141;
t201 = t110 * t141;
t191 = t137 * t141;
t158 = qJD(2) * (t137 * t209 + t191);
t162 = t136 * t148 * t183;
t200 = (t138 * t162 + t146 * t158) / t131 ^ 2;
t199 = t119 * t140;
t197 = t124 * t142;
t196 = t126 * t148;
t194 = t136 * t137;
t147 = t149 ^ 2;
t192 = t136 * t147;
t189 = t141 * t149;
t184 = qJD(1) * t148;
t182 = qJD(2) * t143;
t108 = t110 * t192 + 0.1e1;
t179 = 0.2e1 / t108 ^ 2 * (-t192 * t205 + (t141 * t147 * t182 - t162) * t110);
t178 = 0.2e1 * t205;
t177 = -0.2e1 * t204;
t176 = t124 * t202;
t175 = t110 * t189;
t173 = t129 * t194;
t169 = 0.1e1 + t193;
t168 = t141 * t179;
t167 = t200 * t207;
t166 = t200 * t210;
t163 = t148 * t173;
t161 = t169 * t149;
t160 = t120 * t197 - t199;
t157 = t141 * t181 + t164 * t149;
t123 = -t143 * t185 + t190;
t116 = t169 * t148 * t129;
t114 = 0.1e1 / t117;
t106 = 0.1e1 / t108;
t102 = (t206 * t141 * t126 - t127 * t163) * t149;
t101 = -t143 * t196 + t195 + (t126 * t143 - t127 * t186) * t116;
t99 = -t169 * t166 + (qJD(1) * t161 + t158 * t210) * t129;
t1 = [t137 * t149 * t167 + (qJD(2) * t161 - t184 * t191) * t129, t99, 0, 0, 0, 0; (t109 * t168 + (-t109 * t182 + (qJD(1) * t102 + t98) * t201) * t106) * t148 + (t110 * t168 * t102 + (-((t103 * t163 + t206 * t182 + t167) * t126 + (t166 * t194 - t203 + (t203 + (t207 - t209) * t181) * t129) * t127) * t175 + (-t110 * t182 + t141 * t178) * t102 + (-t109 + ((-t146 + t147) * t127 * t173 + t206 * t174) * t110) * t141 * qJD(1)) * t106) * t149 (t101 * t201 - t109 * t143) * t149 * t179 + ((-t109 * t184 + (-qJD(2) * t101 - t98) * t149 * t110) * t143 + (-t109 * t180 - (-t127 * t148 * t99 - t208 * t126 + (-qJD(2) * t126 + t103 * t196 - t127 * t183) * t116) * t175 + (t110 * t184 + t149 * t178) * t101 - ((t99 - t183) * t126 + ((-t116 * t148 + 0.1e1) * qJD(2) + (t116 - t148) * t103) * t127) * t110 * t188) * t141) * t106, 0, 0, 0, 0; 0.2e1 * (t119 * t159 + t123 * t198) * t204 + (0.2e1 * t123 * t176 - t165 * t119 * t185 + t157 * t199 + (-t165 * t124 * t187 + t123 * t104 + t105 * t159 - t157 * t197) * t120) * t114, t160 * t177 * t189 + (t160 * t143 * t180 + (-t160 * t184 + ((-qJD(5) * t119 - 0.2e1 * t176) * t142 + (-t104 * t142 + (-qJD(5) * t124 + t105) * t140) * t120) * t149) * t141) * t114, 0, 0, t177 + 0.2e1 * (-t104 * t114 * t120 + (-t114 * t202 - t120 * t204) * t124) * t124, 0;];
JaD_rot  = t1;
