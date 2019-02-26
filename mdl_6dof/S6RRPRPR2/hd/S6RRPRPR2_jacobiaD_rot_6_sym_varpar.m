% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:38:16
% EndTime: 2019-02-26 21:38:17
% DurationCPUTime: 0.78s
% Computational Cost: add. (5302->95), mult. (3810->207), div. (753->12), fcn. (4455->9), ass. (0->95)
t177 = sin(qJ(1));
t172 = qJ(2) + pkin(10) + qJ(4);
t170 = sin(t172);
t166 = 0.1e1 / t170 ^ 2;
t171 = cos(t172);
t169 = t171 ^ 2;
t223 = t166 * t169;
t198 = 0.1e1 + t223;
t238 = t177 * t198;
t174 = t177 ^ 2;
t163 = t174 * t223 + 0.1e1;
t161 = 0.1e1 / t163;
t165 = 0.1e1 / t170;
t179 = cos(qJ(1));
t210 = qJD(1) * t179;
t199 = t171 * t210;
t173 = qJD(2) + qJD(4);
t218 = t173 * t177;
t201 = t166 * t218;
t135 = ((t170 * t218 - t199) * t165 + t169 * t201) * t161;
t237 = -t135 + t218;
t194 = qJD(1) * t170 + qJD(6);
t217 = t173 * t179;
t236 = -t171 * t217 + t194 * t177;
t214 = t177 * t171;
t160 = atan2(-t214, t170);
t151 = cos(t160);
t150 = sin(t160);
t203 = t150 * t214;
t145 = t151 * t170 - t203;
t142 = 0.1e1 / t145;
t178 = cos(qJ(6));
t213 = t177 * t178;
t176 = sin(qJ(6));
t215 = t176 * t179;
t157 = t170 * t215 + t213;
t153 = 0.1e1 / t157;
t143 = 0.1e1 / t145 ^ 2;
t154 = 0.1e1 / t157 ^ 2;
t235 = t161 - 0.1e1;
t227 = t151 * t171;
t130 = (-t135 * t177 + t173) * t227 + (t170 * t237 - t199) * t150;
t234 = t130 * t142 * t143;
t195 = qJD(6) * t170 + qJD(1);
t190 = t195 * t179;
t140 = t176 * t190 + t178 * t236;
t212 = t178 * t179;
t216 = t176 * t177;
t156 = -t170 * t212 + t216;
t152 = t156 ^ 2;
t149 = t152 * t154 + 0.1e1;
t226 = t154 * t156;
t141 = -t176 * t236 + t178 * t190;
t231 = t141 * t153 * t154;
t233 = (t140 * t226 - t152 * t231) / t149 ^ 2;
t168 = t171 * t169;
t224 = t165 * t171;
t188 = t173 * (-t165 * t166 * t168 - t224);
t221 = t169 * t177;
t192 = t210 * t221;
t232 = (t166 * t192 + t174 * t188) / t163 ^ 2;
t230 = t143 * t171;
t229 = t143 * t179;
t228 = t150 * t177;
t225 = t156 * t176;
t175 = t179 ^ 2;
t222 = t169 * t175;
t220 = t170 * t173;
t219 = t171 * t173;
t211 = qJD(1) * t177;
t138 = t143 * t222 + 0.1e1;
t209 = 0.2e1 * (-t222 * t234 + (-t170 * t175 * t219 - t192) * t143) / t138 ^ 2;
t208 = 0.2e1 * t234;
t207 = 0.2e1 * t233;
t206 = -0.2e1 * t232;
t205 = t171 * t232;
t204 = t171 * t229;
t202 = t165 * t221;
t197 = t171 * t209;
t196 = 0.2e1 * t156 * t231;
t193 = t151 * t161 * t165 * t169;
t191 = t198 * t179;
t189 = t153 * t178 + t154 * t225;
t187 = t189 * t179;
t159 = -t170 * t216 + t212;
t158 = t170 * t213 + t215;
t147 = 0.1e1 / t149;
t146 = t161 * t238;
t136 = 0.1e1 / t138;
t134 = (t235 * t171 * t150 + t177 * t193) * t179;
t132 = t170 * t228 + t227 + (-t150 * t170 - t151 * t214) * t146;
t131 = t206 * t238 + (qJD(1) * t191 + 0.2e1 * t177 * t188) * t161;
t128 = t171 * t187 * t207 + (t187 * t220 + (t189 * t211 + ((qJD(6) * t153 + t196) * t176 + (-t140 * t176 + (-qJD(6) * t156 + t141) * t178) * t154) * t179) * t171) * t147;
t127 = (t132 * t230 + t142 * t170) * t179 * t209 + ((t142 * t211 + (t132 * t173 + t130) * t229) * t170 + (-t142 * t217 - (-t131 * t151 * t177 + t237 * t150 + (t135 * t228 - t150 * t173 - t151 * t210) * t146) * t204 + (t143 * t211 + t179 * t208) * t132 - ((-t131 + t210) * t150 + ((t146 * t177 - 0.1e1) * t173 + (-t146 + t177) * t135) * t151) * t170 * t229) * t171) * t136;
t1 = [0.2e1 * t165 * t179 * t205 + (t173 * t191 + t211 * t224) * t161, t131, 0, t131, 0, 0; (t142 * t197 + (t142 * t220 + (qJD(1) * t134 + t130) * t230) * t136) * t177 + (t143 * t197 * t134 + (-((-0.2e1 * t205 + t220 + (-t135 * t202 - t220) * t161) * t150 + (t202 * t206 - t135 * t171 + (-t168 * t201 + (t135 - 0.2e1 * t218) * t171) * t161) * t151) * t204 + (t143 * t220 + t171 * t208) * t134 + (-t142 + ((t174 - t175) * t193 + t235 * t203) * t143) * t171 * qJD(1)) * t136) * t179, t127, 0, t127, 0, 0; (-t153 * t158 + t159 * t226) * t207 + (t159 * t196 + (-t159 * t140 - t158 * t141 + t195 * t156 * t213 - (-t173 * t214 - t194 * t179) * t225) * t154 + (t194 * t212 + (-t195 * t176 + t178 * t219) * t177) * t153) * t147, t128, 0, t128, 0, -0.2e1 * t233 + 0.2e1 * (t140 * t147 * t154 + (-t147 * t231 - t154 * t233) * t156) * t156;];
JaD_rot  = t1;
