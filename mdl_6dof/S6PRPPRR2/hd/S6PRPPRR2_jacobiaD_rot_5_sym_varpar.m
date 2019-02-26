% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPPRR2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:13
% EndTime: 2019-02-26 19:45:13
% DurationCPUTime: 0.50s
% Computational Cost: add. (1747->59), mult. (5333->135), div. (281->12), fcn. (6885->13), ass. (0->71)
t181 = sin(pkin(11));
t184 = cos(pkin(11));
t188 = sin(qJ(2));
t190 = cos(qJ(2));
t177 = t188 * t181 - t190 * t184;
t186 = cos(pkin(6));
t198 = t177 * t186;
t201 = t190 * t181 + t188 * t184;
t174 = t201 * t186;
t182 = sin(pkin(10));
t185 = cos(pkin(10));
t159 = t185 * t174 - t182 * t177;
t183 = sin(pkin(6));
t173 = t201 * t183;
t143 = atan2(-t159, t173);
t138 = sin(t143);
t139 = cos(t143);
t132 = -t138 * t159 + t139 * t173;
t129 = 0.1e1 / t132;
t162 = t182 * t198 - t185 * t201;
t187 = sin(qJ(5));
t189 = cos(qJ(5));
t208 = t182 * t183;
t149 = -t162 * t187 + t189 * t208;
t145 = 0.1e1 / t149;
t169 = 0.1e1 / t173;
t130 = 0.1e1 / t132 ^ 2;
t146 = 0.1e1 / t149 ^ 2;
t170 = 0.1e1 / t173 ^ 2;
t156 = t159 ^ 2;
t142 = t156 * t170 + 0.1e1;
t140 = 0.1e1 / t142;
t168 = qJD(2) * t198;
t176 = t201 * qJD(2);
t151 = -t185 * t168 - t182 * t176;
t172 = t177 * t183;
t167 = qJD(2) * t172;
t211 = t159 * t170;
t123 = (-t151 * t169 - t167 * t211) * t140;
t203 = -t138 * t173 - t139 * t159;
t120 = t203 * t123 - t138 * t151 - t139 * t167;
t216 = t120 * t129 * t130;
t148 = t162 * t189 + t187 * t208;
t144 = t148 ^ 2;
t135 = t144 * t146 + 0.1e1;
t175 = t177 * qJD(2);
t197 = t186 * t176;
t152 = t185 * t175 + t182 * t197;
t137 = t149 * qJD(5) + t152 * t189;
t212 = t146 * t148;
t204 = qJD(5) * t148;
t136 = -t152 * t187 - t204;
t213 = t136 * t145 * t146;
t215 = (t137 * t212 - t144 * t213) / t135 ^ 2;
t202 = -t182 * t174 - t185 * t177;
t214 = t130 * t202;
t210 = t159 * t172;
t209 = t167 * t169 * t170;
t153 = t182 * t168 - t185 * t176;
t200 = t145 * t189 + t187 * t212;
t158 = -t182 * t201 - t185 * t198;
t199 = -t158 * t169 - t170 * t210;
t166 = t183 * t176;
t157 = t202 ^ 2;
t150 = t182 * t175 - t185 * t197;
t133 = 0.1e1 / t135;
t127 = t157 * t130 + 0.1e1;
t124 = t199 * t140;
t121 = t203 * t124 - t138 * t158 - t139 * t172;
t119 = -0.2e1 * t199 / t142 ^ 2 * (t151 * t211 + t156 * t209) + (-0.2e1 * t209 * t210 - t150 * t169 + (-t151 * t172 - t158 * t167 - t159 * t166) * t170) * t140;
t1 = [0, t119, 0, 0, 0, 0; 0, 0.2e1 * (t121 * t214 - t129 * t162) / t127 ^ 2 * (t153 * t214 - t157 * t216) + (t152 * t129 + (-t162 * t120 - t121 * t153) * t130 + (0.2e1 * t121 * t216 + (-(-t119 * t159 - t124 * t151 - t166 + (-t124 * t173 - t158) * t123) * t139 - (-t119 * t173 + t124 * t167 - t150 + (t124 * t159 + t172) * t123) * t138) * t130) * t202) / t127, 0, 0, 0, 0; 0, 0.2e1 * t200 * t202 * t215 + (-t200 * t153 + ((qJD(5) * t145 + 0.2e1 * t148 * t213) * t187 + (-t137 * t187 + (t136 - t204) * t189) * t146) * t202) * t133, 0, 0, -0.2e1 * t215 + 0.2e1 * (t133 * t137 * t146 + (-t133 * t213 - t146 * t215) * t148) * t148, 0;];
JaD_rot  = t1;
