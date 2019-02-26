% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:51:42
% EndTime: 2019-02-26 20:51:43
% DurationCPUTime: 0.72s
% Computational Cost: add. (3326->96), mult. (2734->202), div. (498->12), fcn. (3199->9), ass. (0->96)
t164 = pkin(10) + qJ(3);
t161 = sin(t164);
t157 = t161 ^ 2;
t162 = cos(t164);
t159 = 0.1e1 / t162 ^ 2;
t213 = t157 * t159;
t168 = sin(qJ(1));
t231 = 0.2e1 * t168;
t230 = t161 * t213;
t166 = t168 ^ 2;
t151 = t166 * t213 + 0.1e1;
t149 = 0.1e1 / t151;
t158 = 0.1e1 / t162;
t169 = cos(qJ(1));
t203 = qJD(1) * t169;
t191 = t161 * t203;
t201 = qJD(3) * t168;
t124 = (-(-t162 * t201 - t191) * t158 + t201 * t213) * t149;
t229 = t124 - t201;
t163 = pkin(11) + qJ(5) + qJ(6);
t155 = cos(t163);
t205 = t169 * t155;
t154 = sin(t163);
t209 = t168 * t154;
t144 = t162 * t205 + t209;
t207 = t168 * t161;
t148 = atan2(-t207, -t162);
t147 = cos(t148);
t146 = sin(t148);
t194 = t146 * t207;
t134 = -t147 * t162 - t194;
t131 = 0.1e1 / t134;
t138 = 0.1e1 / t144;
t132 = 0.1e1 / t134 ^ 2;
t139 = 0.1e1 / t144 ^ 2;
t228 = -0.2e1 * t161;
t227 = t149 - 0.1e1;
t215 = t147 * t161;
t117 = (-t124 * t168 + qJD(3)) * t215 + (t229 * t162 - t191) * t146;
t226 = t117 * t131 * t132;
t165 = qJD(5) + qJD(6);
t179 = t162 * t209 + t205;
t200 = qJD(3) * t169;
t190 = t161 * t200;
t122 = t179 * qJD(1) - t144 * t165 + t154 * t190;
t206 = t169 * t154;
t208 = t168 * t155;
t143 = t162 * t206 - t208;
t137 = t143 ^ 2;
t130 = t137 * t139 + 0.1e1;
t218 = t139 * t143;
t184 = -qJD(1) * t162 + t165;
t185 = t162 * t165 - qJD(1);
t123 = -t185 * t206 + (t184 * t168 - t190) * t155;
t224 = t123 * t138 * t139;
t225 = (-t122 * t218 - t137 * t224) / t130 ^ 2;
t223 = t124 * t161;
t222 = t132 * t161;
t221 = t132 * t169;
t211 = t158 * t161;
t178 = qJD(3) * (t158 * t230 + t211);
t182 = t157 * t168 * t203;
t220 = (t159 * t182 + t166 * t178) / t151 ^ 2;
t219 = t138 * t154;
t217 = t143 * t155;
t216 = t146 * t168;
t214 = t157 * t158;
t167 = t169 ^ 2;
t212 = t157 * t167;
t210 = t161 * t169;
t204 = qJD(1) * t168;
t202 = qJD(3) * t162;
t127 = t132 * t212 + 0.1e1;
t199 = 0.2e1 * (-t212 * t226 + (t161 * t167 * t202 - t182) * t132) / t127 ^ 2;
t198 = 0.2e1 * t226;
t197 = -0.2e1 * t225;
t196 = t132 * t210;
t195 = t143 * t224;
t193 = t149 * t214;
t189 = 0.1e1 + t213;
t188 = t161 * t199;
t187 = t220 * t228;
t186 = t220 * t231;
t183 = t168 * t193;
t181 = t189 * t169;
t180 = t139 * t217 - t219;
t177 = t161 * t201 + t184 * t169;
t142 = -t162 * t208 + t206;
t136 = t189 * t168 * t149;
t128 = 0.1e1 / t130;
t125 = 0.1e1 / t127;
t121 = (t227 * t161 * t146 - t147 * t183) * t169;
t120 = -t162 * t216 + t215 + (t146 * t162 - t147 * t207) * t136;
t118 = -t189 * t186 + (qJD(1) * t181 + t178 * t231) * t149;
t115 = t197 + 0.2e1 * (-t122 * t139 * t128 + (-t128 * t224 - t139 * t225) * t143) * t143;
t1 = [t169 * t158 * t187 + (qJD(3) * t181 - t204 * t211) * t149, 0, t118, 0, 0, 0; (t131 * t188 + (-t131 * t202 + (qJD(1) * t121 + t117) * t222) * t125) * t168 + (t132 * t188 * t121 + (-((t124 * t183 + t227 * t202 + t187) * t146 + (t186 * t214 - t223 + (t223 + (t228 - t230) * t201) * t149) * t147) * t196 + (-t132 * t202 + t161 * t198) * t121 + (-t131 + ((-t166 + t167) * t147 * t193 + t227 * t194) * t132) * t161 * qJD(1)) * t125) * t169, 0 (t120 * t222 - t131 * t162) * t169 * t199 + ((-t131 * t204 + (-qJD(3) * t120 - t117) * t221) * t162 + (-t131 * t200 - (-t118 * t147 * t168 - t229 * t146 + (-qJD(3) * t146 + t124 * t216 - t147 * t203) * t136) * t196 + (t132 * t204 + t169 * t198) * t120 - ((t118 - t203) * t146 + ((-t136 * t168 + 0.1e1) * qJD(3) + (t136 - t168) * t124) * t147) * t162 * t221) * t161) * t125, 0, 0, 0; 0.2e1 * (t138 * t179 + t142 * t218) * t225 + (0.2e1 * t142 * t195 - t185 * t138 * t208 + t177 * t219 + (-t185 * t143 * t209 + t142 * t122 + t123 * t179 - t177 * t217) * t139) * t128, 0, t180 * t197 * t210 + (t180 * t162 * t200 + (-t180 * t204 + ((-t138 * t165 - 0.2e1 * t195) * t155 + (-t122 * t155 + (-t143 * t165 + t123) * t154) * t139) * t169) * t161) * t128, 0, t115, t115;];
JaD_rot  = t1;
