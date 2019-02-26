% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:50:49
% EndTime: 2019-02-26 20:50:50
% DurationCPUTime: 0.70s
% Computational Cost: add. (2446->95), mult. (2734->206), div. (498->12), fcn. (3199->9), ass. (0->94)
t158 = sin(qJ(3));
t153 = 0.1e1 / t158 ^ 2;
t159 = cos(qJ(3));
t156 = t159 ^ 2;
t197 = t153 * t156;
t151 = qJ(1) + pkin(10);
t146 = sin(t151);
t179 = 0.1e1 + t197;
t219 = t146 * t179;
t218 = t159 * t197;
t157 = qJ(5) + qJ(6);
t148 = sin(t157);
t149 = cos(t157);
t150 = qJD(5) + qJD(6);
t176 = t150 * t158 + qJD(1);
t191 = qJD(3) * t159;
t217 = t176 * t148 - t149 * t191;
t216 = t148 * t191 + t176 * t149;
t201 = t146 * t159;
t139 = atan2(-t201, t158);
t138 = cos(t139);
t137 = sin(t139);
t184 = t137 * t201;
t128 = t138 * t158 - t184;
t124 = 0.1e1 / t128;
t147 = cos(t151);
t200 = t148 * t158;
t134 = t146 * t149 + t147 * t200;
t130 = 0.1e1 / t134;
t152 = 0.1e1 / t158;
t125 = 0.1e1 / t128 ^ 2;
t131 = 0.1e1 / t134 ^ 2;
t144 = t146 ^ 2;
t143 = t144 * t197 + 0.1e1;
t140 = 0.1e1 / t143;
t215 = t140 - 0.1e1;
t175 = qJD(1) * t158 + t150;
t170 = t175 * t149;
t114 = t146 * t170 + t217 * t147;
t199 = t149 * t158;
t133 = t146 * t148 - t147 * t199;
t129 = t133 ^ 2;
t122 = t129 * t131 + 0.1e1;
t205 = t131 * t133;
t171 = t175 * t148;
t115 = -t146 * t171 + t216 * t147;
t212 = t115 * t130 * t131;
t214 = (t114 * t205 - t129 * t212) / t122 ^ 2;
t194 = qJD(1) * t159;
t182 = t147 * t194;
t192 = qJD(3) * t158;
t193 = qJD(3) * t146;
t116 = ((t146 * t192 - t182) * t152 + t193 * t197) * t140;
t203 = t138 * t159;
t110 = (-t116 * t146 + qJD(3)) * t203 + (-t182 + (-t116 + t193) * t158) * t137;
t213 = t110 * t124 * t125;
t211 = t116 * t137;
t210 = t116 * t159;
t168 = qJD(3) * (-t159 - t218) * t152;
t195 = qJD(1) * t147;
t173 = t146 * t156 * t195;
t209 = (t144 * t168 + t153 * t173) / t143 ^ 2;
t208 = t125 * t147;
t207 = t125 * t159;
t127 = t140 * t219;
t206 = t127 * t146;
t204 = t137 * t158;
t145 = t147 ^ 2;
t202 = t145 * t156;
t198 = t152 * t156;
t196 = qJD(1) * t146;
t119 = t125 * t202 + 0.1e1;
t190 = 0.2e1 * (-t202 * t213 + (-t145 * t158 * t191 - t173) * t125) / t119 ^ 2;
t189 = 0.2e1 * t214;
t188 = 0.2e1 * t213;
t187 = -0.2e1 * t209;
t186 = t159 * t209;
t185 = t147 * t207;
t183 = t146 * t198;
t178 = t159 * t190;
t177 = 0.2e1 * t133 * t212;
t174 = t140 * t183;
t172 = t179 * t147;
t169 = t130 * t149 + t148 * t205;
t167 = t169 * t159;
t136 = -t146 * t200 + t147 * t149;
t135 = t146 * t199 + t147 * t148;
t120 = 0.1e1 / t122;
t117 = 0.1e1 / t119;
t113 = (t215 * t159 * t137 + t138 * t174) * t147;
t112 = t146 * t204 + t203 + (-t138 * t201 - t204) * t127;
t111 = t187 * t219 + (qJD(1) * t172 + 0.2e1 * t146 * t168) * t140;
t107 = -0.2e1 * t214 + 0.2e1 * (t114 * t120 * t131 + (-t120 * t212 - t131 * t214) * t133) * t133;
t1 = [0.2e1 * t147 * t152 * t186 + (t146 * t152 * t194 + qJD(3) * t172) * t140, 0, t111, 0, 0, 0; (t124 * t178 + (t124 * t192 + (qJD(1) * t113 + t110) * t207) * t117) * t146 + (t125 * t178 * t113 + (-((-t116 * t174 - t215 * t192 - 0.2e1 * t186) * t137 + (t183 * t187 - t210 + (t210 + (-0.2e1 * t159 - t218) * t193) * t140) * t138) * t185 + (t125 * t192 + t159 * t188) * t113 + (-t124 + ((t144 - t145) * t140 * t138 * t198 + t215 * t184) * t125) * t194) * t117) * t147, 0 (t112 * t207 + t124 * t158) * t147 * t190 + ((t124 * t196 + (qJD(3) * t112 + t110) * t208) * t158 + (-t147 * qJD(3) * t124 - (-t111 * t138 * t146 + t137 * t193 + t206 * t211 - t211 + (-qJD(3) * t137 - t138 * t195) * t127) * t185 + (t125 * t196 + t147 * t188) * t112 - ((-t111 + t195) * t137 + ((-0.1e1 + t206) * qJD(3) + (-t127 + t146) * t116) * t138) * t158 * t208) * t159) * t117, 0, 0, 0; (-t130 * t135 + t136 * t205) * t189 + (t136 * t177 + (-t136 * t114 - t135 * t115 + (t216 * t146 + t147 * t171) * t133) * t131 + (-t217 * t146 + t147 * t170) * t130) * t120, 0, t147 * t167 * t189 + (t167 * t196 + (t169 * t192 + ((t130 * t150 + t177) * t148 + (-t114 * t148 + (-t133 * t150 + t115) * t149) * t131) * t159) * t147) * t120, 0, t107, t107;];
JaD_rot  = t1;
