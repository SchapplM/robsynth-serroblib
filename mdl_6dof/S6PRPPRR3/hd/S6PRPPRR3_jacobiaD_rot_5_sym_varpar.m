% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPPRR3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:50
% EndTime: 2019-02-26 19:45:51
% DurationCPUTime: 0.51s
% Computational Cost: add. (1747->63), mult. (5333->142), div. (281->12), fcn. (6885->13), ass. (0->74)
t183 = sin(pkin(11));
t184 = sin(pkin(10));
t188 = sin(qJ(2));
t190 = cos(qJ(2));
t220 = cos(pkin(10));
t221 = cos(pkin(6));
t205 = t221 * t220;
t198 = -t184 * t190 - t188 * t205;
t222 = t198 * t183;
t185 = sin(pkin(6));
t186 = cos(pkin(11));
t175 = (t183 * t188 + t186 * t190) * t185;
t209 = t184 * t221;
t181 = -t188 * t209 + t220 * t190;
t200 = t184 * t188 - t190 * t205;
t159 = -t200 * t186 - t222;
t144 = atan2(-t159, t175);
t139 = sin(t144);
t140 = cos(t144);
t133 = -t139 * t159 + t140 * t175;
t130 = 0.1e1 / t133;
t199 = -t220 * t188 - t190 * t209;
t164 = t181 * t186 - t183 * t199;
t187 = sin(qJ(5));
t189 = cos(qJ(5));
t211 = t184 * t185;
t150 = t164 * t189 - t187 * t211;
t146 = 0.1e1 / t150;
t172 = 0.1e1 / t175;
t131 = 0.1e1 / t133 ^ 2;
t147 = 0.1e1 / t150 ^ 2;
t173 = 0.1e1 / t175 ^ 2;
t156 = t159 ^ 2;
t143 = t156 * t173 + 0.1e1;
t141 = 0.1e1 / t143;
t177 = t200 * qJD(2);
t197 = t198 * t186;
t151 = qJD(2) * t197 - t177 * t183;
t176 = (t183 * t190 - t186 * t188) * t185;
t168 = qJD(2) * t176;
t214 = t159 * t173;
t124 = (-t151 * t172 + t168 * t214) * t141;
t204 = -t139 * t175 - t140 * t159;
t121 = t204 * t124 - t139 * t151 + t140 * t168;
t219 = t121 * t130 * t131;
t149 = t164 * t187 + t189 * t211;
t145 = t149 ^ 2;
t136 = t145 * t147 + 0.1e1;
t178 = t199 * qJD(2);
t179 = t181 * qJD(2);
t154 = t178 * t186 + t179 * t183;
t137 = t150 * qJD(5) + t154 * t187;
t215 = t147 * t149;
t210 = qJD(5) * t149;
t138 = t154 * t189 - t210;
t216 = t138 * t146 * t147;
t218 = (t137 * t215 - t145 * t216) / t136 ^ 2;
t207 = t181 * t183 + t186 * t199;
t217 = t131 * t207;
t213 = t159 * t176;
t212 = t168 * t172 * t173;
t153 = t178 * t183 - t179 * t186;
t202 = -t146 * t187 + t189 * t215;
t158 = -t200 * t183 + t197;
t201 = -t158 * t172 + t173 * t213;
t167 = qJD(2) * t175;
t157 = t207 ^ 2;
t152 = qJD(2) * t222 + t177 * t186;
t134 = 0.1e1 / t136;
t128 = t157 * t131 + 0.1e1;
t125 = t201 * t141;
t122 = t204 * t125 - t139 * t158 + t140 * t176;
t120 = -0.2e1 * t201 / t143 ^ 2 * (t151 * t214 - t156 * t212) + (-0.2e1 * t212 * t213 - t152 * t172 + (t151 * t176 + t158 * t168 - t159 * t167) * t173) * t141;
t1 = [0, t120, 0, 0, 0, 0; 0, 0.2e1 * (t122 * t217 + t130 * t164) / t128 ^ 2 * (t153 * t217 - t157 * t219) + (-t154 * t130 + (t121 * t164 - t122 * t153) * t131 + (0.2e1 * t122 * t219 + (-(-t120 * t159 - t125 * t151 - t167 + (-t125 * t175 - t158) * t124) * t140 - (-t120 * t175 - t125 * t168 - t152 + (t125 * t159 - t176) * t124) * t139) * t131) * t207) / t128, 0, 0, 0, 0; 0, 0.2e1 * t202 * t207 * t218 + (-t202 * t153 + ((qJD(5) * t146 + 0.2e1 * t149 * t216) * t189 + (-t137 * t189 + (-t138 + t210) * t187) * t147) * t207) * t134, 0, 0, -0.2e1 * t218 + 0.2e1 * (t134 * t137 * t147 + (-t134 * t216 - t147 * t218) * t149) * t149, 0;];
JaD_rot  = t1;
