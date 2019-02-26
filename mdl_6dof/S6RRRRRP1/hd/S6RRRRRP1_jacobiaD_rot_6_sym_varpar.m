% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP1
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
% Datum: 2019-02-26 22:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRP1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:39:34
% EndTime: 2019-02-26 22:39:34
% DurationCPUTime: 0.76s
% Computational Cost: add. (7770->97), mult. (5101->208), div. (1026->12), fcn. (5942->9), ass. (0->95)
t177 = sin(qJ(1));
t239 = 0.2e1 * t177;
t174 = t177 ^ 2;
t173 = qJ(2) + qJ(3) + qJ(4);
t170 = sin(t173);
t165 = t170 ^ 2;
t171 = cos(t173);
t167 = 0.1e1 / t171 ^ 2;
t224 = t165 * t167;
t161 = t174 * t224 + 0.1e1;
t159 = 0.1e1 / t161;
t166 = 0.1e1 / t171;
t179 = cos(qJ(1));
t210 = qJD(1) * t179;
t200 = t170 * t210;
t172 = qJD(2) + qJD(3) + qJD(4);
t218 = t172 * t177;
t203 = t167 * t218;
t133 = (-(-t171 * t218 - t200) * t166 + t165 * t203) * t159;
t238 = t133 - t218;
t178 = cos(qJ(5));
t212 = t179 * t178;
t176 = sin(qJ(5));
t215 = t177 * t176;
t157 = t171 * t212 + t215;
t216 = t177 * t170;
t158 = atan2(-t216, -t171);
t149 = cos(t158);
t148 = sin(t158);
t205 = t148 * t216;
t143 = -t149 * t171 - t205;
t140 = 0.1e1 / t143;
t151 = 0.1e1 / t157;
t141 = 0.1e1 / t143 ^ 2;
t152 = 0.1e1 / t157 ^ 2;
t237 = t159 - 0.1e1;
t229 = t149 * t170;
t128 = (-t133 * t177 + t172) * t229 + (t238 * t171 - t200) * t148;
t236 = t128 * t140 * t141;
t188 = t171 * t215 + t212;
t217 = t172 * t179;
t201 = t170 * t217;
t138 = t188 * qJD(1) - t157 * qJD(5) + t176 * t201;
t213 = t179 * t176;
t214 = t177 * t178;
t156 = t171 * t213 - t214;
t150 = t156 ^ 2;
t147 = t150 * t152 + 0.1e1;
t227 = t152 * t156;
t194 = -qJD(1) * t171 + qJD(5);
t195 = qJD(5) * t171 - qJD(1);
t139 = -t195 * t213 + (t194 * t177 - t201) * t178;
t233 = t139 * t151 * t152;
t235 = (-t138 * t227 - t150 * t233) / t147 ^ 2;
t164 = t170 * t165;
t221 = t166 * t170;
t187 = t172 * (t164 * t166 * t167 + t221);
t222 = t165 * t177;
t192 = t210 * t222;
t234 = (t167 * t192 + t174 * t187) / t161 ^ 2;
t232 = t141 * t170;
t231 = t141 * t179;
t230 = t148 * t177;
t228 = t151 * t176;
t226 = t156 * t178;
t225 = t165 * t166;
t175 = t179 ^ 2;
t223 = t165 * t175;
t220 = t170 * t179;
t219 = t171 * t172;
t211 = qJD(1) * t177;
t136 = t141 * t223 + 0.1e1;
t209 = 0.2e1 * (-t223 * t236 + (t170 * t175 * t219 - t192) * t141) / t136 ^ 2;
t208 = 0.2e1 * t236;
t207 = -0.2e1 * t235;
t206 = t141 * t220;
t204 = t156 * t233;
t199 = 0.1e1 + t224;
t198 = t170 * t209;
t197 = -0.2e1 * t170 * t234;
t196 = t234 * t239;
t193 = t149 * t159 * t225;
t191 = t179 * t199;
t190 = t194 * t179;
t189 = t152 * t226 - t228;
t155 = -t171 * t214 + t213;
t145 = 0.1e1 / t147;
t144 = t199 * t177 * t159;
t134 = 0.1e1 / t136;
t132 = (t237 * t170 * t148 - t177 * t193) * t179;
t130 = -t171 * t230 + t229 + (t148 * t171 - t149 * t216) * t144;
t129 = -t199 * t196 + (qJD(1) * t191 + t187 * t239) * t159;
t126 = t189 * t207 * t220 + (t189 * t171 * t217 + (-t189 * t211 + ((-qJD(5) * t151 - 0.2e1 * t204) * t178 + (-t138 * t178 + (-qJD(5) * t156 + t139) * t176) * t152) * t179) * t170) * t145;
t125 = (t130 * t232 - t140 * t171) * t179 * t209 + ((-t140 * t211 + (-t130 * t172 - t128) * t231) * t171 + (-t140 * t217 - (-t129 * t149 * t177 - t238 * t148 + (t133 * t230 - t148 * t172 - t149 * t210) * t144) * t206 + (t141 * t211 + t179 * t208) * t130 - ((t129 - t210) * t148 + ((-t144 * t177 + 0.1e1) * t172 + (t144 - t177) * t133) * t149) * t171 * t231) * t170) * t134;
t1 = [t179 * t166 * t197 + (t172 * t191 - t211 * t221) * t159, t129, t129, t129, 0, 0; (t140 * t198 + (-t140 * t219 + (qJD(1) * t132 + t128) * t232) * t134) * t177 + (t141 * t198 * t132 + (-((t197 - t219 + (t133 * t166 * t222 + t219) * t159) * t148 + (t196 * t225 - t133 * t170 + (-t164 * t203 + (t133 - 0.2e1 * t218) * t170) * t159) * t149) * t206 + (-t141 * t219 + t170 * t208) * t132 + (-t140 + ((-t174 + t175) * t193 + t237 * t205) * t141) * t170 * qJD(1)) * t134) * t179, t125, t125, t125, 0, 0; 0.2e1 * (t151 * t188 + t155 * t227) * t235 + (0.2e1 * t155 * t204 - t195 * t151 * t214 + (t172 * t216 + t190) * t228 + (t155 * t138 + t188 * t139 - t190 * t226 - (t170 * t172 * t178 + t195 * t176) * t156 * t177) * t152) * t145, t126, t126, t126, t207 + 0.2e1 * (-t138 * t152 * t145 + (-t145 * t233 - t152 * t235) * t156) * t156, 0;];
JaD_rot  = t1;
