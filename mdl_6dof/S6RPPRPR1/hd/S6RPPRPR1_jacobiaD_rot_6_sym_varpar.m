% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:25:30
% EndTime: 2019-02-26 20:25:30
% DurationCPUTime: 0.76s
% Computational Cost: add. (3592->98), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->94)
t150 = pkin(10) + qJ(4);
t144 = sin(t150);
t137 = t144 ^ 2;
t147 = cos(t150);
t140 = 0.1e1 / t147 ^ 2;
t196 = t137 * t140;
t151 = qJ(1) + pkin(9);
t145 = sin(t151);
t213 = 0.2e1 * t145;
t212 = t144 * t196;
t149 = pkin(11) + qJ(6);
t143 = sin(t149);
t146 = cos(t149);
t148 = cos(t151);
t187 = t147 * t148;
t126 = t143 * t145 + t146 * t187;
t190 = t145 * t144;
t129 = atan2(-t190, -t147);
t128 = cos(t129);
t127 = sin(t129);
t176 = t127 * t190;
t113 = -t128 * t147 - t176;
t110 = 0.1e1 / t113;
t120 = 0.1e1 / t126;
t139 = 0.1e1 / t147;
t111 = 0.1e1 / t113 ^ 2;
t121 = 0.1e1 / t126 ^ 2;
t211 = -0.2e1 * t144;
t138 = t145 ^ 2;
t132 = t138 * t196 + 0.1e1;
t130 = 0.1e1 / t132;
t210 = t130 - 0.1e1;
t185 = qJD(1) * t148;
t173 = t144 * t185;
t183 = qJD(4) * t147;
t184 = qJD(4) * t145;
t104 = (-(-t145 * t183 - t173) * t139 + t184 * t196) * t130;
t198 = t128 * t144;
t99 = (-t104 * t145 + qJD(4)) * t198 + (-t173 + (t104 - t184) * t147) * t127;
t209 = t110 * t111 * t99;
t188 = t145 * t147;
t160 = t143 * t188 + t146 * t148;
t182 = qJD(4) * t148;
t172 = t144 * t182;
t105 = t160 * qJD(1) - t126 * qJD(6) + t143 * t172;
t189 = t145 * t146;
t125 = t143 * t187 - t189;
t119 = t125 ^ 2;
t118 = t119 * t121 + 0.1e1;
t200 = t121 * t125;
t166 = -qJD(1) * t147 + qJD(6);
t167 = qJD(6) * t147 - qJD(1);
t192 = t143 * t148;
t106 = -t167 * t192 + (t145 * t166 - t172) * t146;
t205 = t106 * t120 * t121;
t208 = (-t105 * t200 - t119 * t205) / t118 ^ 2;
t207 = t104 * t127;
t206 = t104 * t144;
t204 = t111 * t144;
t194 = t139 * t144;
t159 = qJD(4) * (t139 * t212 + t194);
t164 = t137 * t145 * t185;
t203 = (t138 * t159 + t140 * t164) / t132 ^ 2;
t171 = 0.1e1 + t196;
t117 = t171 * t145 * t130;
t202 = t117 * t145;
t201 = t120 * t143;
t199 = t125 * t146;
t197 = t137 * t139;
t142 = t148 ^ 2;
t195 = t137 * t142;
t191 = t144 * t148;
t186 = qJD(1) * t145;
t109 = t111 * t195 + 0.1e1;
t181 = 0.2e1 / t109 ^ 2 * (-t195 * t209 + (t142 * t144 * t183 - t164) * t111);
t180 = 0.2e1 * t209;
t179 = -0.2e1 * t208;
t178 = t125 * t205;
t177 = t111 * t191;
t175 = t130 * t197;
t170 = t144 * t181;
t169 = t203 * t211;
t168 = t203 * t213;
t165 = t145 * t175;
t163 = t171 * t148;
t162 = t166 * t148;
t161 = t121 * t199 - t201;
t124 = -t146 * t188 + t192;
t115 = 0.1e1 / t118;
t107 = 0.1e1 / t109;
t103 = (t127 * t144 * t210 - t128 * t165) * t148;
t102 = -t127 * t188 + t198 + (t127 * t147 - t128 * t190) * t117;
t100 = -t171 * t168 + (qJD(1) * t163 + t159 * t213) * t130;
t1 = [t139 * t148 * t169 + (qJD(4) * t163 - t186 * t194) * t130, 0, 0, t100, 0, 0; (t110 * t170 + (-t110 * t183 + (qJD(1) * t103 + t99) * t204) * t107) * t145 + (t111 * t170 * t103 + (-((t104 * t165 + t183 * t210 + t169) * t127 + (t168 * t197 - t206 + (t206 + (t211 - t212) * t184) * t130) * t128) * t177 + (-t111 * t183 + t144 * t180) * t103 + (-t110 + ((-t138 + t142) * t128 * t175 + t210 * t176) * t111) * t144 * qJD(1)) * t107) * t148, 0, 0 (t102 * t204 - t110 * t147) * t148 * t181 + ((-t110 * t186 + (-qJD(4) * t102 - t99) * t148 * t111) * t147 + (-t110 * t182 - (-t100 * t128 * t145 + t127 * t184 + t202 * t207 - t207 + (-qJD(4) * t127 - t128 * t185) * t117) * t177 + (t111 * t186 + t148 * t180) * t102 - ((t100 - t185) * t127 + ((0.1e1 - t202) * qJD(4) + (t117 - t145) * t104) * t128) * t111 * t187) * t144) * t107, 0, 0; 0.2e1 * (t120 * t160 + t124 * t200) * t208 + (0.2e1 * t124 * t178 - t167 * t120 * t189 + (t144 * t184 + t162) * t201 + (t124 * t105 + t160 * t106 - t162 * t199 - (qJD(4) * t144 * t146 + t143 * t167) * t125 * t145) * t121) * t115, 0, 0, t161 * t179 * t191 + (t161 * t147 * t182 + (-t161 * t186 + ((-qJD(6) * t120 - 0.2e1 * t178) * t146 + (-t105 * t146 + (-qJD(6) * t125 + t106) * t143) * t121) * t148) * t144) * t115, 0, t179 + 0.2e1 * (-t105 * t115 * t121 + (-t115 * t205 - t121 * t208) * t125) * t125;];
JaD_rot  = t1;
