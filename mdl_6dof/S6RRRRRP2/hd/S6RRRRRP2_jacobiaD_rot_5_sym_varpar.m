% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRP2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:40:10
% EndTime: 2019-02-26 22:40:11
% DurationCPUTime: 0.77s
% Computational Cost: add. (7770->97), mult. (5101->208), div. (1026->12), fcn. (5942->9), ass. (0->95)
t179 = sin(qJ(1));
t241 = 0.2e1 * t179;
t176 = t179 ^ 2;
t175 = qJ(2) + qJ(3) + qJ(4);
t172 = sin(t175);
t167 = t172 ^ 2;
t173 = cos(t175);
t169 = 0.1e1 / t173 ^ 2;
t226 = t167 * t169;
t163 = t176 * t226 + 0.1e1;
t161 = 0.1e1 / t163;
t168 = 0.1e1 / t173;
t181 = cos(qJ(1));
t212 = qJD(1) * t181;
t202 = t172 * t212;
t174 = qJD(2) + qJD(3) + qJD(4);
t220 = t174 * t179;
t205 = t169 * t220;
t135 = (-(-t173 * t220 - t202) * t168 + t167 * t205) * t161;
t240 = t135 - t220;
t180 = cos(qJ(5));
t214 = t180 * t181;
t178 = sin(qJ(5));
t216 = t179 * t178;
t159 = t173 * t214 + t216;
t217 = t179 * t172;
t160 = atan2(-t217, -t173);
t151 = cos(t160);
t150 = sin(t160);
t206 = t150 * t217;
t145 = -t151 * t173 - t206;
t142 = 0.1e1 / t145;
t153 = 0.1e1 / t159;
t143 = 0.1e1 / t145 ^ 2;
t154 = 0.1e1 / t159 ^ 2;
t239 = t161 - 0.1e1;
t231 = t151 * t172;
t130 = (-t135 * t179 + t174) * t231 + (t240 * t173 - t202) * t150;
t238 = t130 * t142 * t143;
t166 = t172 * t167;
t223 = t168 * t172;
t189 = t174 * (t166 * t168 * t169 + t223);
t224 = t167 * t179;
t194 = t212 * t224;
t237 = (t169 * t194 + t176 * t189) / t163 ^ 2;
t196 = -qJD(1) * t173 + qJD(5);
t197 = qJD(5) * t173 - qJD(1);
t219 = t174 * t181;
t203 = t172 * t219;
t218 = t178 * t181;
t141 = -t197 * t218 + (t196 * t179 - t203) * t180;
t236 = t141 * t153 * t154;
t235 = t143 * t172;
t234 = t143 * t181;
t190 = t173 * t216 + t214;
t140 = t190 * qJD(1) - t159 * qJD(5) + t178 * t203;
t215 = t179 * t180;
t158 = t173 * t218 - t215;
t152 = t158 ^ 2;
t149 = t152 * t154 + 0.1e1;
t229 = t154 * t158;
t233 = 0.1e1 / t149 ^ 2 * (-t140 * t229 - t152 * t236);
t232 = t150 * t179;
t230 = t153 * t178;
t228 = t158 * t180;
t227 = t167 * t168;
t177 = t181 ^ 2;
t225 = t167 * t177;
t222 = t172 * t181;
t221 = t173 * t174;
t213 = qJD(1) * t179;
t138 = t143 * t225 + 0.1e1;
t211 = 0.2e1 * (-t225 * t238 + (t172 * t177 * t221 - t194) * t143) / t138 ^ 2;
t210 = 0.2e1 * t238;
t209 = -0.2e1 * t233;
t208 = t158 * t236;
t207 = t143 * t222;
t201 = 0.1e1 + t226;
t200 = t172 * t211;
t199 = -0.2e1 * t172 * t237;
t198 = t237 * t241;
t195 = t151 * t161 * t227;
t193 = t201 * t181;
t192 = t196 * t181;
t191 = t154 * t228 - t230;
t157 = -t173 * t215 + t218;
t147 = 0.1e1 / t149;
t146 = t201 * t179 * t161;
t136 = 0.1e1 / t138;
t134 = (t239 * t172 * t150 - t179 * t195) * t181;
t132 = -t173 * t232 + t231 + (t150 * t173 - t151 * t217) * t146;
t131 = -t201 * t198 + (qJD(1) * t193 + t189 * t241) * t161;
t128 = t191 * t209 * t222 + (t191 * t173 * t219 + (-t191 * t213 + ((-qJD(5) * t153 - 0.2e1 * t208) * t180 + (-t140 * t180 + (-qJD(5) * t158 + t141) * t178) * t154) * t181) * t172) * t147;
t127 = (t132 * t235 - t142 * t173) * t181 * t211 + ((-t142 * t213 + (-t132 * t174 - t130) * t234) * t173 + (-t142 * t219 - (-t131 * t151 * t179 - t240 * t150 + (t135 * t232 - t150 * t174 - t151 * t212) * t146) * t207 + (t143 * t213 + t181 * t210) * t132 - ((t131 - t212) * t150 + ((-t146 * t179 + 0.1e1) * t174 + (t146 - t179) * t135) * t151) * t173 * t234) * t172) * t136;
t1 = [t181 * t168 * t199 + (t174 * t193 - t213 * t223) * t161, t131, t131, t131, 0, 0; (t142 * t200 + (-t142 * t221 + (qJD(1) * t134 + t130) * t235) * t136) * t179 + (t143 * t200 * t134 + (-((t199 - t221 + (t135 * t168 * t224 + t221) * t161) * t150 + (t198 * t227 - t135 * t172 + (-t166 * t205 + (t135 - 0.2e1 * t220) * t172) * t161) * t151) * t207 + (-t143 * t221 + t172 * t210) * t134 + (-t142 + ((-t176 + t177) * t195 + t239 * t206) * t143) * t172 * qJD(1)) * t136) * t181, t127, t127, t127, 0, 0; 0.2e1 * (t153 * t190 + t157 * t229) * t233 + (0.2e1 * t157 * t208 - t197 * t153 * t215 + (t174 * t217 + t192) * t230 + (t157 * t140 + t190 * t141 - t192 * t228 - (t172 * t174 * t180 + t197 * t178) * t158 * t179) * t154) * t147, t128, t128, t128, t209 + 0.2e1 * (-t140 * t154 * t147 + (-t147 * t236 - t154 * t233) * t158) * t158, 0;];
JaD_rot  = t1;
