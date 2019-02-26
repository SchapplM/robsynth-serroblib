% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:04:06
% EndTime: 2019-02-26 20:04:06
% DurationCPUTime: 0.40s
% Computational Cost: add. (1419->57), mult. (2577->133), div. (441->14), fcn. (3313->11), ass. (0->69)
t177 = sin(qJ(2));
t178 = cos(qJ(2));
t175 = sin(pkin(11));
t208 = cos(pkin(6));
t193 = t175 * t208;
t207 = cos(pkin(11));
t162 = -t177 * t193 + t207 * t178;
t189 = t208 * t207;
t158 = t175 * t177 - t178 * t189;
t176 = sin(pkin(6));
t197 = t176 * t178;
t148 = atan2(-t158, -t197);
t146 = sin(t148);
t147 = cos(t148);
t133 = -t146 * t158 - t147 * t197;
t130 = 0.1e1 / t133;
t168 = qJ(3) + pkin(12) + qJ(5);
t166 = sin(t168);
t167 = cos(t168);
t198 = t175 * t176;
t145 = t162 * t167 + t166 * t198;
t141 = 0.1e1 / t145;
t172 = 0.1e1 / t178;
t131 = 0.1e1 / t133 ^ 2;
t142 = 0.1e1 / t145 ^ 2;
t173 = 0.1e1 / t178 ^ 2;
t160 = t175 * t178 + t177 * t189;
t153 = t160 * qJD(2);
t199 = t173 * t177;
t194 = t158 * t199;
t156 = t158 ^ 2;
t170 = 0.1e1 / t176 ^ 2;
t151 = t156 * t170 * t173 + 0.1e1;
t149 = 0.1e1 / t151;
t169 = 0.1e1 / t176;
t201 = t149 * t169;
t125 = (qJD(2) * t194 + t153 * t172) * t201;
t187 = t146 * t197 - t147 * t158;
t195 = t147 * t176 * t177;
t122 = qJD(2) * t195 + t125 * t187 - t146 * t153;
t206 = t122 * t130 * t131;
t144 = t162 * t166 - t167 * t198;
t140 = t144 ^ 2;
t136 = t140 * t142 + 0.1e1;
t185 = -t177 * t207 - t178 * t193;
t154 = t185 * qJD(2);
t171 = qJD(3) + qJD(5);
t190 = t171 * t198 + t154;
t200 = t162 * t171;
t137 = t166 * t190 + t167 * t200;
t202 = t142 * t144;
t138 = -t166 * t200 + t167 * t190;
t203 = t138 * t141 * t142;
t205 = (t137 * t202 - t140 * t203) / t136 ^ 2;
t204 = t131 * t185;
t196 = -0.2e1 * t205;
t188 = -t141 * t166 + t167 * t202;
t186 = t160 * t172 + t194;
t174 = t172 * t173;
t157 = t185 ^ 2;
t155 = t162 * qJD(2);
t152 = qJD(2) * t158;
t134 = 0.1e1 / t136;
t129 = t157 * t131 + 0.1e1;
t126 = t186 * t201;
t123 = t126 * t187 - t146 * t160 + t195;
t121 = (-0.2e1 * t186 / t151 ^ 2 * (qJD(2) * t156 * t174 * t177 + t153 * t158 * t173) * t170 + (t153 * t199 - t152 * t172 + (t160 * t199 + (0.2e1 * t174 * t177 ^ 2 + t172) * t158) * qJD(2)) * t149) * t169;
t119 = t196 + 0.2e1 * (t134 * t137 * t142 + (-t134 * t203 - t142 * t205) * t144) * t144;
t1 = [0, t121, 0, 0, 0, 0; 0, 0.2e1 * (-t123 * t204 - t130 * t162) / t129 ^ 2 * (-t155 * t204 - t157 * t206) + (t154 * t130 + (-t162 * t122 - t123 * t155) * t131 - (0.2e1 * t123 * t206 + (-(qJD(2) * t197 - t121 * t158 - t126 * t153 + (t126 * t197 - t160) * t125) * t147 - (t125 * t126 * t158 + t152 + (t121 * t178 + (-qJD(2) * t126 - t125) * t177) * t176) * t146) * t131) * t185) / t129, 0, 0, 0, 0; 0, -t188 * t185 * t196 + (t188 * t155 - ((-t141 * t171 - 0.2e1 * t144 * t203) * t167 + (t137 * t167 + (-t144 * t171 + t138) * t166) * t142) * t185) * t134, t119, 0, t119, 0;];
JaD_rot  = t1;
