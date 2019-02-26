% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:03:56
% EndTime: 2019-02-26 22:03:57
% DurationCPUTime: 0.74s
% Computational Cost: add. (5302->95), mult. (3810->207), div. (753->12), fcn. (4455->9), ass. (0->95)
t179 = sin(qJ(1));
t174 = qJ(2) + qJ(3) + pkin(10);
t172 = sin(t174);
t168 = 0.1e1 / t172 ^ 2;
t173 = cos(t174);
t171 = t173 ^ 2;
t225 = t168 * t171;
t200 = 0.1e1 + t225;
t240 = t179 * t200;
t176 = t179 ^ 2;
t165 = t176 * t225 + 0.1e1;
t163 = 0.1e1 / t165;
t167 = 0.1e1 / t172;
t181 = cos(qJ(1));
t212 = qJD(1) * t181;
t201 = t173 * t212;
t175 = qJD(2) + qJD(3);
t220 = t175 * t179;
t203 = t168 * t220;
t137 = ((t172 * t220 - t201) * t167 + t171 * t203) * t163;
t239 = -t137 + t220;
t196 = qJD(1) * t172 + qJD(6);
t219 = t175 * t181;
t238 = -t173 * t219 + t196 * t179;
t216 = t179 * t173;
t162 = atan2(-t216, t172);
t153 = cos(t162);
t152 = sin(t162);
t205 = t152 * t216;
t147 = t153 * t172 - t205;
t144 = 0.1e1 / t147;
t180 = cos(qJ(6));
t215 = t179 * t180;
t178 = sin(qJ(6));
t217 = t178 * t181;
t159 = t172 * t217 + t215;
t155 = 0.1e1 / t159;
t145 = 0.1e1 / t147 ^ 2;
t156 = 0.1e1 / t159 ^ 2;
t237 = t163 - 0.1e1;
t229 = t153 * t173;
t132 = (-t137 * t179 + t175) * t229 + (t172 * t239 - t201) * t152;
t236 = t132 * t144 * t145;
t197 = qJD(6) * t172 + qJD(1);
t192 = t197 * t181;
t142 = t178 * t192 + t180 * t238;
t214 = t180 * t181;
t218 = t178 * t179;
t158 = -t172 * t214 + t218;
t154 = t158 ^ 2;
t151 = t154 * t156 + 0.1e1;
t228 = t156 * t158;
t143 = -t178 * t238 + t180 * t192;
t233 = t143 * t155 * t156;
t235 = (t142 * t228 - t154 * t233) / t151 ^ 2;
t170 = t173 * t171;
t226 = t167 * t173;
t190 = t175 * (-t167 * t168 * t170 - t226);
t223 = t171 * t179;
t194 = t212 * t223;
t234 = (t168 * t194 + t176 * t190) / t165 ^ 2;
t232 = t145 * t173;
t231 = t145 * t181;
t230 = t152 * t179;
t227 = t158 * t178;
t177 = t181 ^ 2;
t224 = t171 * t177;
t222 = t172 * t175;
t221 = t173 * t175;
t213 = qJD(1) * t179;
t140 = t145 * t224 + 0.1e1;
t211 = 0.2e1 * (-t224 * t236 + (-t172 * t177 * t221 - t194) * t145) / t140 ^ 2;
t210 = 0.2e1 * t236;
t209 = 0.2e1 * t235;
t208 = -0.2e1 * t234;
t207 = t173 * t234;
t206 = t173 * t231;
t204 = t167 * t223;
t199 = t173 * t211;
t198 = 0.2e1 * t158 * t233;
t195 = t153 * t163 * t167 * t171;
t193 = t200 * t181;
t191 = t155 * t180 + t156 * t227;
t189 = t191 * t181;
t161 = -t172 * t218 + t214;
t160 = t172 * t215 + t217;
t149 = 0.1e1 / t151;
t148 = t163 * t240;
t138 = 0.1e1 / t140;
t136 = (t237 * t173 * t152 + t179 * t195) * t181;
t134 = t172 * t230 + t229 + (-t152 * t172 - t153 * t216) * t148;
t133 = t208 * t240 + (qJD(1) * t193 + 0.2e1 * t179 * t190) * t163;
t130 = t173 * t189 * t209 + (t189 * t222 + (t191 * t213 + ((qJD(6) * t155 + t198) * t178 + (-t142 * t178 + (-qJD(6) * t158 + t143) * t180) * t156) * t181) * t173) * t149;
t129 = (t134 * t232 + t144 * t172) * t181 * t211 + ((t144 * t213 + (t134 * t175 + t132) * t231) * t172 + (-t144 * t219 - (-t133 * t153 * t179 + t239 * t152 + (t137 * t230 - t152 * t175 - t153 * t212) * t148) * t206 + (t145 * t213 + t181 * t210) * t134 - ((-t133 + t212) * t152 + ((t148 * t179 - 0.1e1) * t175 + (-t148 + t179) * t137) * t153) * t172 * t231) * t173) * t138;
t1 = [0.2e1 * t167 * t181 * t207 + (t175 * t193 + t213 * t226) * t163, t133, t133, 0, 0, 0; (t144 * t199 + (t144 * t222 + (qJD(1) * t136 + t132) * t232) * t138) * t179 + (t145 * t199 * t136 + (-((-0.2e1 * t207 + t222 + (-t137 * t204 - t222) * t163) * t152 + (t204 * t208 - t137 * t173 + (-t170 * t203 + (t137 - 0.2e1 * t220) * t173) * t163) * t153) * t206 + (t145 * t222 + t173 * t210) * t136 + (-t144 + ((t176 - t177) * t195 + t237 * t205) * t145) * t173 * qJD(1)) * t138) * t181, t129, t129, 0, 0, 0; (-t155 * t160 + t161 * t228) * t209 + (t161 * t198 + (-t161 * t142 - t160 * t143 + t197 * t158 * t215 - (-t175 * t216 - t196 * t181) * t227) * t156 + (t196 * t214 + (-t197 * t178 + t180 * t221) * t179) * t155) * t149, t130, t130, 0, 0, -0.2e1 * t235 + 0.2e1 * (t142 * t149 * t156 + (-t149 * t233 - t156 * t235) * t158) * t158;];
JaD_rot  = t1;
