% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:50
% EndTime: 2019-02-26 21:37:51
% DurationCPUTime: 0.78s
% Computational Cost: add. (6073->97), mult. (3810->205), div. (753->12), fcn. (4455->9), ass. (0->96)
t179 = sin(qJ(1));
t240 = 0.2e1 * t179;
t177 = t179 ^ 2;
t174 = qJ(2) + pkin(10) + qJ(4);
t170 = sin(t174);
t166 = t170 ^ 2;
t171 = cos(t174);
t168 = 0.1e1 / t171 ^ 2;
t225 = t166 * t168;
t161 = t177 * t225 + 0.1e1;
t159 = 0.1e1 / t161;
t167 = 0.1e1 / t171;
t180 = cos(qJ(1));
t211 = qJD(1) * t180;
t201 = t170 * t211;
t176 = qJD(2) + qJD(4);
t219 = t176 * t179;
t204 = t168 * t219;
t133 = (-(-t171 * t219 - t201) * t167 + t166 * t204) * t159;
t239 = t133 - t219;
t175 = pkin(11) + qJ(6);
t173 = cos(t175);
t213 = t180 * t173;
t172 = sin(t175);
t216 = t179 * t172;
t155 = t171 * t213 + t216;
t217 = t179 * t170;
t158 = atan2(-t217, -t171);
t157 = cos(t158);
t156 = sin(t158);
t205 = t156 * t217;
t143 = -t157 * t171 - t205;
t140 = 0.1e1 / t143;
t149 = 0.1e1 / t155;
t141 = 0.1e1 / t143 ^ 2;
t150 = 0.1e1 / t155 ^ 2;
t238 = t159 - 0.1e1;
t227 = t157 * t170;
t128 = (-t133 * t179 + t176) * t227 + (t239 * t171 - t201) * t156;
t237 = t128 * t140 * t141;
t190 = t171 * t216 + t213;
t218 = t176 * t180;
t202 = t170 * t218;
t134 = t190 * qJD(1) - t155 * qJD(6) + t172 * t202;
t214 = t180 * t172;
t215 = t179 * t173;
t154 = t171 * t214 - t215;
t148 = t154 ^ 2;
t147 = t148 * t150 + 0.1e1;
t230 = t150 * t154;
t195 = -qJD(1) * t171 + qJD(6);
t196 = qJD(6) * t171 - qJD(1);
t135 = -t196 * t214 + (t195 * t179 - t202) * t173;
t235 = t135 * t149 * t150;
t236 = (-t134 * t230 - t148 * t235) / t147 ^ 2;
t165 = t170 * t166;
t222 = t167 * t170;
t189 = t176 * (t165 * t167 * t168 + t222);
t223 = t166 * t179;
t193 = t211 * t223;
t234 = (t168 * t193 + t177 * t189) / t161 ^ 2;
t233 = t141 * t170;
t232 = t141 * t180;
t231 = t149 * t172;
t229 = t154 * t173;
t228 = t156 * t179;
t226 = t166 * t167;
t178 = t180 ^ 2;
t224 = t166 * t178;
t221 = t170 * t180;
t220 = t171 * t176;
t212 = qJD(1) * t179;
t138 = t141 * t224 + 0.1e1;
t210 = 0.2e1 * (-t224 * t237 + (t170 * t178 * t220 - t193) * t141) / t138 ^ 2;
t209 = 0.2e1 * t237;
t208 = -0.2e1 * t236;
t207 = t141 * t221;
t206 = t154 * t235;
t200 = 0.1e1 + t225;
t199 = t170 * t210;
t198 = -0.2e1 * t170 * t234;
t197 = t234 * t240;
t194 = t157 * t159 * t226;
t192 = t200 * t180;
t191 = t150 * t229 - t231;
t188 = t176 * t217 + t195 * t180;
t153 = -t171 * t215 + t214;
t145 = 0.1e1 / t147;
t144 = t200 * t179 * t159;
t136 = 0.1e1 / t138;
t132 = (t238 * t170 * t156 - t179 * t194) * t180;
t131 = -t171 * t228 + t227 + (t156 * t171 - t157 * t217) * t144;
t129 = -t200 * t197 + (qJD(1) * t192 + t189 * t240) * t159;
t126 = t191 * t208 * t221 + (t191 * t171 * t218 + (-t191 * t212 + ((-qJD(6) * t149 - 0.2e1 * t206) * t173 + (-t134 * t173 + (-qJD(6) * t154 + t135) * t172) * t150) * t180) * t170) * t145;
t125 = (t131 * t233 - t140 * t171) * t180 * t210 + ((-t140 * t212 + (-t131 * t176 - t128) * t232) * t171 + (-t140 * t218 - (-t129 * t157 * t179 - t239 * t156 + (t133 * t228 - t156 * t176 - t157 * t211) * t144) * t207 + (t141 * t212 + t180 * t209) * t131 - ((t129 - t211) * t156 + ((-t144 * t179 + 0.1e1) * t176 + (t144 - t179) * t133) * t157) * t171 * t232) * t170) * t136;
t1 = [t180 * t167 * t198 + (t176 * t192 - t212 * t222) * t159, t129, 0, t129, 0, 0; (t140 * t199 + (-t140 * t220 + (qJD(1) * t132 + t128) * t233) * t136) * t179 + (t141 * t199 * t132 + (-((t198 - t220 + (t133 * t167 * t223 + t220) * t159) * t156 + (t197 * t226 - t133 * t170 + (-t165 * t204 + (t133 - 0.2e1 * t219) * t170) * t159) * t157) * t207 + (-t141 * t220 + t170 * t209) * t132 + (-t140 + ((-t177 + t178) * t194 + t238 * t205) * t141) * t170 * qJD(1)) * t136) * t180, t125, 0, t125, 0, 0; 0.2e1 * (t149 * t190 + t153 * t230) * t236 + (0.2e1 * t153 * t206 - t196 * t149 * t215 + t188 * t231 + (-t196 * t154 * t216 + t153 * t134 + t135 * t190 - t188 * t229) * t150) * t145, t126, 0, t126, 0, t208 + 0.2e1 * (-t134 * t150 * t145 + (-t145 * t235 - t150 * t236) * t154) * t154;];
JaD_rot  = t1;
