% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRP1_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:13
% EndTime: 2019-02-26 20:15:13
% DurationCPUTime: 0.46s
% Computational Cost: add. (1157->57), mult. (2577->133), div. (441->14), fcn. (3313->11), ass. (0->69)
t172 = sin(qJ(2));
t173 = cos(qJ(2));
t170 = sin(pkin(11));
t203 = cos(pkin(6));
t188 = t170 * t203;
t202 = cos(pkin(11));
t157 = -t172 * t188 + t202 * t173;
t184 = t203 * t202;
t153 = t170 * t172 - t173 * t184;
t171 = sin(pkin(6));
t192 = t171 * t173;
t143 = atan2(-t153, -t192);
t141 = sin(t143);
t142 = cos(t143);
t128 = -t141 * t153 - t142 * t192;
t125 = 0.1e1 / t128;
t169 = qJ(3) + qJ(4);
t161 = sin(t169);
t162 = cos(t169);
t193 = t170 * t171;
t140 = t157 * t162 + t161 * t193;
t136 = 0.1e1 / t140;
t166 = 0.1e1 / t173;
t126 = 0.1e1 / t128 ^ 2;
t137 = 0.1e1 / t140 ^ 2;
t167 = 0.1e1 / t173 ^ 2;
t155 = t170 * t173 + t172 * t184;
t148 = t155 * qJD(2);
t194 = t167 * t172;
t189 = t153 * t194;
t151 = t153 ^ 2;
t164 = 0.1e1 / t171 ^ 2;
t146 = t151 * t164 * t167 + 0.1e1;
t144 = 0.1e1 / t146;
t163 = 0.1e1 / t171;
t196 = t144 * t163;
t120 = (qJD(2) * t189 + t148 * t166) * t196;
t182 = t141 * t192 - t142 * t153;
t190 = t142 * t171 * t172;
t117 = qJD(2) * t190 + t182 * t120 - t141 * t148;
t201 = t117 * t125 * t126;
t180 = -t202 * t172 - t173 * t188;
t200 = t126 * t180;
t139 = t157 * t161 - t162 * t193;
t135 = t139 ^ 2;
t132 = t135 * t137 + 0.1e1;
t149 = t180 * qJD(2);
t165 = qJD(3) + qJD(4);
t185 = t165 * t193 + t149;
t195 = t157 * t165;
t133 = t185 * t161 + t162 * t195;
t197 = t137 * t139;
t134 = -t161 * t195 + t185 * t162;
t198 = t134 * t136 * t137;
t199 = 0.1e1 / t132 ^ 2 * (t133 * t197 - t135 * t198);
t191 = -0.2e1 * t199;
t183 = -t136 * t161 + t162 * t197;
t181 = t155 * t166 + t189;
t168 = t166 * t167;
t152 = t180 ^ 2;
t150 = t157 * qJD(2);
t147 = qJD(2) * t153;
t129 = 0.1e1 / t132;
t124 = t126 * t152 + 0.1e1;
t121 = t181 * t196;
t118 = t182 * t121 - t141 * t155 + t190;
t116 = (-0.2e1 * t181 / t146 ^ 2 * (qJD(2) * t151 * t168 * t172 + t148 * t153 * t167) * t164 + (t148 * t194 - t147 * t166 + (t155 * t194 + (0.2e1 * t168 * t172 ^ 2 + t166) * t153) * qJD(2)) * t144) * t163;
t114 = t191 + 0.2e1 * (t129 * t133 * t137 + (-t129 * t198 - t137 * t199) * t139) * t139;
t1 = [0, t116, 0, 0, 0, 0; 0, 0.2e1 * (-t118 * t200 - t125 * t157) / t124 ^ 2 * (-t150 * t200 - t152 * t201) + (t149 * t125 + (-t157 * t117 - t118 * t150) * t126 - (0.2e1 * t118 * t201 + (-(qJD(2) * t192 - t116 * t153 - t121 * t148 + (t121 * t192 - t155) * t120) * t142 - (t120 * t121 * t153 + t147 + (t116 * t173 + (-qJD(2) * t121 - t120) * t172) * t171) * t141) * t126) * t180) / t124, 0, 0, 0, 0; 0, -t183 * t180 * t191 + (t183 * t150 - ((-t136 * t165 - 0.2e1 * t139 * t198) * t162 + (t133 * t162 + (-t139 * t165 + t134) * t161) * t137) * t180) * t129, t114, t114, 0, 0;];
JaD_rot  = t1;
