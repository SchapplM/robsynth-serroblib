% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:10
% EndTime: 2019-02-26 21:55:11
% DurationCPUTime: 0.74s
% Computational Cost: add. (3714->96), mult. (2949->202), div. (516->12), fcn. (3430->9), ass. (0->96)
t172 = qJ(2) + pkin(11);
t168 = sin(t172);
t164 = t168 ^ 2;
t169 = cos(t172);
t166 = 0.1e1 / t169 ^ 2;
t220 = t164 * t166;
t175 = sin(qJ(1));
t238 = 0.2e1 * t175;
t237 = t168 * t220;
t173 = t175 ^ 2;
t158 = t173 * t220 + 0.1e1;
t156 = 0.1e1 / t158;
t165 = 0.1e1 / t169;
t176 = cos(qJ(1));
t210 = qJD(1) * t176;
t198 = t168 * t210;
t208 = qJD(2) * t175;
t131 = (-(-t169 * t208 - t198) * t165 + t208 * t220) * t156;
t236 = t131 - t208;
t171 = qJ(4) + qJ(5) + qJ(6);
t162 = cos(t171);
t212 = t176 * t162;
t161 = sin(t171);
t216 = t175 * t161;
t151 = t169 * t212 + t216;
t214 = t175 * t168;
t155 = atan2(-t214, -t169);
t154 = cos(t155);
t153 = sin(t155);
t201 = t153 * t214;
t140 = -t154 * t169 - t201;
t135 = 0.1e1 / t140;
t145 = 0.1e1 / t151;
t136 = 0.1e1 / t140 ^ 2;
t146 = 0.1e1 / t151 ^ 2;
t235 = -0.2e1 * t168;
t234 = t156 - 0.1e1;
t222 = t154 * t168;
t124 = (-t131 * t175 + qJD(2)) * t222 + (t236 * t169 - t198) * t153;
t233 = t124 * t135 * t136;
t170 = qJD(4) + qJD(5) + qJD(6);
t186 = t169 * t216 + t212;
t207 = qJD(2) * t176;
t197 = t168 * t207;
t129 = qJD(1) * t186 - t151 * t170 + t161 * t197;
t213 = t176 * t161;
t215 = t175 * t162;
t150 = t169 * t213 - t215;
t144 = t150 ^ 2;
t141 = t144 * t146 + 0.1e1;
t225 = t146 * t150;
t191 = -qJD(1) * t169 + t170;
t192 = t169 * t170 - qJD(1);
t130 = -t192 * t213 + (t175 * t191 - t197) * t162;
t231 = t130 * t145 * t146;
t232 = (-t129 * t225 - t144 * t231) / t141 ^ 2;
t230 = t131 * t168;
t229 = t136 * t168;
t228 = t136 * t176;
t218 = t165 * t168;
t185 = qJD(2) * (t165 * t237 + t218);
t189 = t164 * t175 * t210;
t227 = (t166 * t189 + t173 * t185) / t158 ^ 2;
t226 = t145 * t161;
t224 = t150 * t162;
t223 = t153 * t175;
t221 = t164 * t165;
t174 = t176 ^ 2;
t219 = t164 * t174;
t217 = t168 * t176;
t211 = qJD(1) * t175;
t209 = qJD(2) * t169;
t134 = t136 * t219 + 0.1e1;
t206 = 0.2e1 * (-t219 * t233 + (t168 * t174 * t209 - t189) * t136) / t134 ^ 2;
t205 = 0.2e1 * t233;
t204 = -0.2e1 * t232;
t203 = t136 * t217;
t202 = t150 * t231;
t200 = t156 * t221;
t196 = 0.1e1 + t220;
t195 = t168 * t206;
t194 = t227 * t235;
t193 = t227 * t238;
t190 = t175 * t200;
t188 = t196 * t176;
t187 = t224 * t146 - t226;
t184 = t168 * t208 + t176 * t191;
t149 = -t169 * t215 + t213;
t143 = t196 * t175 * t156;
t138 = 0.1e1 / t141;
t132 = 0.1e1 / t134;
t128 = (t153 * t168 * t234 - t154 * t190) * t176;
t127 = -t169 * t223 + t222 + (t153 * t169 - t154 * t214) * t143;
t125 = -t196 * t193 + (qJD(1) * t188 + t185 * t238) * t156;
t122 = t204 + 0.2e1 * (-t129 * t146 * t138 + (-t138 * t231 - t146 * t232) * t150) * t150;
t1 = [t176 * t165 * t194 + (qJD(2) * t188 - t211 * t218) * t156, t125, 0, 0, 0, 0; (t135 * t195 + (-t135 * t209 + (qJD(1) * t128 + t124) * t229) * t132) * t175 + (t136 * t195 * t128 + (-((t131 * t190 + t209 * t234 + t194) * t153 + (t193 * t221 - t230 + (t230 + (t235 - t237) * t208) * t156) * t154) * t203 + (-t136 * t209 + t168 * t205) * t128 + (-t135 + ((-t173 + t174) * t154 * t200 + t234 * t201) * t136) * t168 * qJD(1)) * t132) * t176 (t127 * t229 - t135 * t169) * t176 * t206 + ((-t135 * t211 + (-qJD(2) * t127 - t124) * t228) * t169 + (-t135 * t207 - (-t125 * t154 * t175 - t236 * t153 + (-qJD(2) * t153 + t131 * t223 - t154 * t210) * t143) * t203 + (t136 * t211 + t176 * t205) * t127 - ((t125 - t210) * t153 + ((-t143 * t175 + 0.1e1) * qJD(2) + (t143 - t175) * t131) * t154) * t169 * t228) * t168) * t132, 0, 0, 0, 0; 0.2e1 * (t145 * t186 + t149 * t225) * t232 + (0.2e1 * t149 * t202 - t192 * t145 * t215 + t184 * t226 + (-t150 * t192 * t216 + t149 * t129 + t130 * t186 - t184 * t224) * t146) * t138, t187 * t204 * t217 + (t187 * t169 * t207 + (-t187 * t211 + ((-t145 * t170 - 0.2e1 * t202) * t162 + (-t129 * t162 + (-t150 * t170 + t130) * t161) * t146) * t176) * t168) * t138, 0, t122, t122, t122;];
JaD_rot  = t1;
