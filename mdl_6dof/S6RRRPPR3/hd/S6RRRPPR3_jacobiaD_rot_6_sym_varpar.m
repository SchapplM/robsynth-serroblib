% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:04:32
% EndTime: 2019-02-26 22:04:33
% DurationCPUTime: 0.73s
% Computational Cost: add. (3360->97), mult. (3810->208), div. (753->12), fcn. (4455->9), ass. (0->97)
t177 = sin(qJ(1));
t175 = qJ(2) + qJ(3);
t170 = sin(t175);
t166 = 0.1e1 / t170 ^ 2;
t171 = cos(t175);
t169 = t171 ^ 2;
t224 = t166 * t169;
t198 = 0.1e1 + t224;
t239 = t177 * t198;
t173 = t177 ^ 2;
t162 = t173 * t224 + 0.1e1;
t159 = 0.1e1 / t162;
t165 = 0.1e1 / t170;
t179 = cos(qJ(1));
t211 = qJD(1) * t179;
t199 = t171 * t211;
t172 = qJD(2) + qJD(3);
t219 = t172 * t177;
t202 = t166 * t219;
t133 = ((t170 * t219 - t199) * t165 + t169 * t202) * t159;
t238 = -t133 + t219;
t216 = t177 * t171;
t158 = atan2(-t216, t170);
t157 = cos(t158);
t156 = sin(t158);
t204 = t156 * t216;
t143 = t157 * t170 - t204;
t140 = 0.1e1 / t143;
t178 = cos(qJ(6));
t213 = t178 * t179;
t201 = t170 * t213;
t176 = sin(qJ(6));
t215 = t177 * t176;
t155 = t201 - t215;
t149 = 0.1e1 / t155;
t141 = 0.1e1 / t143 ^ 2;
t150 = 0.1e1 / t155 ^ 2;
t237 = t159 - 0.1e1;
t226 = t157 * t171;
t128 = (-t133 * t177 + t172) * t226 + (t170 * t238 - t199) * t156;
t236 = t128 * t140 * t141;
t194 = qJD(1) * t170 + qJD(6);
t218 = t172 * t179;
t200 = t171 * t218;
t137 = -qJD(6) * t201 - t176 * t200 - t178 * t211 + t194 * t215;
t214 = t177 * t178;
t217 = t176 * t179;
t154 = t170 * t217 + t214;
t148 = t154 ^ 2;
t147 = t148 * t150 + 0.1e1;
t229 = t150 * t154;
t195 = qJD(6) * t170 + qJD(1);
t138 = -t195 * t217 + (-t194 * t177 + t200) * t178;
t234 = t138 * t149 * t150;
t235 = (-t137 * t229 - t148 * t234) / t147 ^ 2;
t168 = t171 * t169;
t225 = t165 * t171;
t188 = t172 * (-t165 * t166 * t168 - t225);
t222 = t169 * t177;
t192 = t211 * t222;
t233 = (t166 * t192 + t173 * t188) / t162 ^ 2;
t232 = t141 * t171;
t231 = t141 * t179;
t230 = t149 * t176;
t228 = t154 * t178;
t227 = t156 * t177;
t174 = t179 ^ 2;
t223 = t169 * t174;
t221 = t170 * t172;
t220 = t171 * t172;
t212 = qJD(1) * t177;
t136 = t141 * t223 + 0.1e1;
t210 = 0.2e1 * (-t223 * t236 + (-t170 * t174 * t220 - t192) * t141) / t136 ^ 2;
t209 = 0.2e1 * t236;
t208 = 0.2e1 * t235;
t207 = -0.2e1 * t233;
t206 = t171 * t233;
t205 = t171 * t231;
t203 = t165 * t222;
t197 = t171 * t210;
t196 = 0.2e1 * t154 * t234;
t193 = t157 * t159 * t165 * t169;
t191 = t198 * t179;
t190 = t194 * t179;
t189 = t150 * t228 - t230;
t187 = t189 * t179;
t153 = -t170 * t214 - t217;
t152 = -t170 * t215 + t213;
t145 = 0.1e1 / t147;
t144 = t159 * t239;
t134 = 0.1e1 / t136;
t132 = (t237 * t171 * t156 + t177 * t193) * t179;
t130 = t170 * t227 + t226 + (-t156 * t170 - t157 * t216) * t144;
t129 = t207 * t239 + (qJD(1) * t191 + 0.2e1 * t177 * t188) * t159;
t126 = t171 * t187 * t208 + (t187 * t221 + (t189 * t212 + ((qJD(6) * t149 + t196) * t178 + (t137 * t178 + (qJD(6) * t154 - t138) * t176) * t150) * t179) * t171) * t145;
t125 = (t130 * t232 + t140 * t170) * t179 * t210 + ((t140 * t212 + (t130 * t172 + t128) * t231) * t170 + (-t140 * t218 - (-t129 * t157 * t177 + t238 * t156 + (t133 * t227 - t156 * t172 - t157 * t211) * t144) * t205 + (t141 * t212 + t179 * t209) * t130 - ((-t129 + t211) * t156 + ((t144 * t177 - 0.1e1) * t172 + (-t144 + t177) * t133) * t157) * t170 * t231) * t171) * t134;
t1 = [0.2e1 * t165 * t179 * t206 + (t172 * t191 + t212 * t225) * t159, t129, t129, 0, 0, 0; (t140 * t197 + (t140 * t221 + (qJD(1) * t132 + t128) * t232) * t134) * t177 + (t141 * t197 * t132 + (-((-0.2e1 * t206 + t221 + (-t133 * t203 - t221) * t159) * t156 + (t203 * t207 - t133 * t171 + (-t168 * t202 + (t133 - 0.2e1 * t219) * t171) * t159) * t157) * t205 + (t141 * t221 + t171 * t209) * t132 + (-t140 + ((t173 - t174) * t193 + t237 * t204) * t141) * t171 * qJD(1)) * t134) * t179, t125, t125, 0, 0, 0; (-t149 * t152 + t153 * t229) * t208 + (t153 * t196 - t195 * t149 * t214 + (-t172 * t216 - t190) * t230 + (t153 * t137 - t152 * t138 + t190 * t228 - (t195 * t176 - t178 * t220) * t154 * t177) * t150) * t145, t126, t126, 0, 0, -0.2e1 * t235 + 0.2e1 * (-t137 * t145 * t150 + (-t145 * t234 - t150 * t235) * t154) * t154;];
JaD_rot  = t1;
