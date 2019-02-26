% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPP1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:52
% EndTime: 2019-02-26 21:34:52
% DurationCPUTime: 0.74s
% Computational Cost: add. (2618->95), mult. (2519->203), div. (480->12), fcn. (2968->9), ass. (0->93)
t147 = qJ(2) + pkin(9);
t143 = sin(t147);
t138 = t143 ^ 2;
t145 = cos(t147);
t140 = 0.1e1 / t145 ^ 2;
t195 = t138 * t140;
t150 = sin(qJ(1));
t212 = 0.2e1 * t150;
t211 = t143 * t195;
t148 = t150 ^ 2;
t133 = t148 * t195 + 0.1e1;
t131 = 0.1e1 / t133;
t139 = 0.1e1 / t145;
t151 = cos(qJ(1));
t185 = qJD(1) * t151;
t173 = t143 * t185;
t183 = qJD(2) * t150;
t105 = (-(-t145 * t183 - t173) * t139 + t183 * t195) * t131;
t210 = t105 - t183;
t146 = qJ(4) + pkin(10);
t144 = cos(t146);
t142 = sin(t146);
t189 = t150 * t142;
t190 = t145 * t151;
t127 = t144 * t190 + t189;
t188 = t150 * t143;
t130 = atan2(-t188, -t145);
t129 = cos(t130);
t128 = sin(t130);
t176 = t128 * t188;
t114 = -t129 * t145 - t176;
t111 = 0.1e1 / t114;
t121 = 0.1e1 / t127;
t112 = 0.1e1 / t114 ^ 2;
t122 = 0.1e1 / t127 ^ 2;
t209 = -0.2e1 * t143;
t208 = t131 - 0.1e1;
t197 = t129 * t143;
t100 = (-t105 * t150 + qJD(2)) * t197 + (t145 * t210 - t173) * t128;
t207 = t100 * t111 * t112;
t161 = t144 * t151 + t145 * t189;
t182 = qJD(2) * t151;
t172 = t143 * t182;
t106 = t161 * qJD(1) - qJD(4) * t127 + t142 * t172;
t187 = t150 * t144;
t126 = t142 * t190 - t187;
t120 = t126 ^ 2;
t119 = t120 * t122 + 0.1e1;
t200 = t122 * t126;
t166 = -qJD(1) * t145 + qJD(4);
t167 = qJD(4) * t145 - qJD(1);
t192 = t142 * t151;
t107 = -t167 * t192 + (t166 * t150 - t172) * t144;
t204 = t107 * t121 * t122;
t206 = (-t106 * t200 - t120 * t204) / t119 ^ 2;
t205 = t105 * t143;
t203 = t112 * t143;
t193 = t139 * t143;
t160 = qJD(2) * (t139 * t211 + t193);
t164 = t138 * t150 * t185;
t202 = (t140 * t164 + t148 * t160) / t133 ^ 2;
t201 = t121 * t142;
t199 = t126 * t144;
t198 = t128 * t150;
t196 = t138 * t139;
t149 = t151 ^ 2;
t194 = t138 * t149;
t191 = t143 * t151;
t186 = qJD(1) * t150;
t184 = qJD(2) * t145;
t110 = t112 * t194 + 0.1e1;
t181 = 0.2e1 / t110 ^ 2 * (-t194 * t207 + (t143 * t149 * t184 - t164) * t112);
t180 = 0.2e1 * t207;
t179 = -0.2e1 * t206;
t178 = t126 * t204;
t177 = t112 * t191;
t175 = t131 * t196;
t171 = 0.1e1 + t195;
t170 = t143 * t181;
t169 = t202 * t209;
t168 = t202 * t212;
t165 = t150 * t175;
t163 = t171 * t151;
t162 = t122 * t199 - t201;
t159 = t143 * t183 + t166 * t151;
t125 = -t145 * t187 + t192;
t118 = t171 * t150 * t131;
t116 = 0.1e1 / t119;
t108 = 0.1e1 / t110;
t104 = (t208 * t143 * t128 - t129 * t165) * t151;
t103 = -t145 * t198 + t197 + (t128 * t145 - t129 * t188) * t118;
t101 = -t171 * t168 + (qJD(1) * t163 + t160 * t212) * t131;
t1 = [t139 * t151 * t169 + (qJD(2) * t163 - t186 * t193) * t131, t101, 0, 0, 0, 0; (t111 * t170 + (-t111 * t184 + (qJD(1) * t104 + t100) * t203) * t108) * t150 + (t112 * t170 * t104 + (-((t105 * t165 + t208 * t184 + t169) * t128 + (t168 * t196 - t205 + (t205 + (t209 - t211) * t183) * t131) * t129) * t177 + (-t112 * t184 + t143 * t180) * t104 + (-t111 + ((-t148 + t149) * t129 * t175 + t208 * t176) * t112) * t143 * qJD(1)) * t108) * t151 (t103 * t203 - t111 * t145) * t151 * t181 + ((-t111 * t186 + (-qJD(2) * t103 - t100) * t151 * t112) * t145 + (-t111 * t182 - (-t101 * t129 * t150 - t210 * t128 + (-qJD(2) * t128 + t105 * t198 - t129 * t185) * t118) * t177 + (t112 * t186 + t151 * t180) * t103 - ((t101 - t185) * t128 + ((-t118 * t150 + 0.1e1) * qJD(2) + (t118 - t150) * t105) * t129) * t112 * t190) * t143) * t108, 0, 0, 0, 0; 0.2e1 * (t121 * t161 + t125 * t200) * t206 + (0.2e1 * t125 * t178 - t167 * t121 * t187 + t159 * t201 + (-t167 * t126 * t189 + t125 * t106 + t107 * t161 - t159 * t199) * t122) * t116, t162 * t179 * t191 + (t162 * t145 * t182 + (-t162 * t186 + ((-qJD(4) * t121 - 0.2e1 * t178) * t144 + (-t106 * t144 + (-qJD(4) * t126 + t107) * t142) * t122) * t151) * t143) * t116, 0, t179 + 0.2e1 * (-t106 * t116 * t122 + (-t116 * t204 - t122 * t206) * t126) * t126, 0, 0;];
JaD_rot  = t1;
