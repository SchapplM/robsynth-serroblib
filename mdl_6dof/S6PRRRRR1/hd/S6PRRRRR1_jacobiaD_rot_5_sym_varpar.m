% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRR1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:36
% EndTime: 2019-02-26 20:18:37
% DurationCPUTime: 0.41s
% Computational Cost: add. (1740->57), mult. (2883->133), div. (459->14), fcn. (3670->11), ass. (0->69)
t176 = sin(qJ(2));
t177 = cos(qJ(2));
t174 = sin(pkin(12));
t207 = cos(pkin(6));
t192 = t174 * t207;
t206 = cos(pkin(12));
t161 = -t176 * t192 + t206 * t177;
t188 = t207 * t206;
t157 = t174 * t176 - t177 * t188;
t175 = sin(pkin(6));
t196 = t175 * t177;
t147 = atan2(-t157, -t196);
t145 = sin(t147);
t146 = cos(t147);
t132 = -t145 * t157 - t196 * t146;
t129 = 0.1e1 / t132;
t168 = qJ(3) + qJ(4) + qJ(5);
t165 = sin(t168);
t166 = cos(t168);
t197 = t174 * t175;
t144 = t161 * t166 + t165 * t197;
t140 = 0.1e1 / t144;
t171 = 0.1e1 / t177;
t130 = 0.1e1 / t132 ^ 2;
t141 = 0.1e1 / t144 ^ 2;
t172 = 0.1e1 / t177 ^ 2;
t159 = t174 * t177 + t176 * t188;
t152 = t159 * qJD(2);
t198 = t172 * t176;
t193 = t157 * t198;
t155 = t157 ^ 2;
t170 = 0.1e1 / t175 ^ 2;
t150 = t155 * t170 * t172 + 0.1e1;
t148 = 0.1e1 / t150;
t169 = 0.1e1 / t175;
t200 = t148 * t169;
t124 = (qJD(2) * t193 + t152 * t171) * t200;
t186 = t145 * t196 - t146 * t157;
t194 = t146 * t175 * t176;
t121 = qJD(2) * t194 + t124 * t186 - t145 * t152;
t205 = t121 * t129 * t130;
t143 = t161 * t165 - t166 * t197;
t139 = t143 ^ 2;
t135 = t139 * t141 + 0.1e1;
t184 = -t176 * t206 - t177 * t192;
t153 = t184 * qJD(2);
t167 = qJD(3) + qJD(4) + qJD(5);
t189 = t167 * t197 + t153;
t199 = t161 * t167;
t136 = t165 * t189 + t166 * t199;
t201 = t141 * t143;
t137 = -t165 * t199 + t166 * t189;
t202 = t137 * t140 * t141;
t204 = (t136 * t201 - t139 * t202) / t135 ^ 2;
t203 = t130 * t184;
t195 = -0.2e1 * t204;
t187 = -t140 * t165 + t166 * t201;
t185 = t159 * t171 + t193;
t173 = t171 * t172;
t156 = t184 ^ 2;
t154 = t161 * qJD(2);
t151 = qJD(2) * t157;
t133 = 0.1e1 / t135;
t128 = t156 * t130 + 0.1e1;
t125 = t185 * t200;
t122 = t125 * t186 - t145 * t159 + t194;
t120 = (-0.2e1 * t185 / t150 ^ 2 * (qJD(2) * t155 * t173 * t176 + t152 * t157 * t172) * t170 + (t152 * t198 - t151 * t171 + (t159 * t198 + (0.2e1 * t173 * t176 ^ 2 + t171) * t157) * qJD(2)) * t148) * t169;
t118 = t195 + 0.2e1 * (t133 * t136 * t141 + (-t133 * t202 - t141 * t204) * t143) * t143;
t1 = [0, t120, 0, 0, 0, 0; 0, 0.2e1 * (-t122 * t203 - t129 * t161) / t128 ^ 2 * (-t154 * t203 - t156 * t205) + (t153 * t129 + (-t161 * t121 - t122 * t154) * t130 - (0.2e1 * t122 * t205 + (-(qJD(2) * t196 - t120 * t157 - t125 * t152 + (t125 * t196 - t159) * t124) * t146 - (t124 * t125 * t157 + t151 + (t120 * t177 + (-qJD(2) * t125 - t124) * t176) * t175) * t145) * t130) * t184) / t128, 0, 0, 0, 0; 0, -t187 * t184 * t195 + (t187 * t154 - ((-t140 * t167 - 0.2e1 * t143 * t202) * t166 + (t136 * t166 + (-t143 * t167 + t137) * t165) * t141) * t184) * t133, t118, t118, t118, 0;];
JaD_rot  = t1;
