% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP10_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:04
% EndTime: 2019-02-26 21:13:06
% DurationCPUTime: 1.11s
% Computational Cost: add. (6985->120), mult. (8378->262), div. (1558->15), fcn. (10537->9), ass. (0->115)
t188 = sin(qJ(1));
t189 = cos(qJ(3));
t236 = t188 * t189;
t186 = qJ(4) + qJ(5);
t178 = sin(t186);
t187 = sin(qJ(3));
t259 = cos(qJ(1));
t219 = t259 * t187;
t179 = cos(t186);
t237 = t188 * t179;
t167 = t178 * t219 + t237;
t235 = t189 * t178;
t155 = atan2(t167, t235);
t152 = cos(t155);
t151 = sin(t155);
t250 = t151 * t167;
t146 = t152 * t235 + t250;
t143 = 0.1e1 / t146;
t166 = t259 * t178 + t187 * t237;
t161 = 0.1e1 / t166;
t175 = 0.1e1 / t178;
t183 = 0.1e1 / t189;
t144 = 0.1e1 / t146 ^ 2;
t162 = 0.1e1 / t166 ^ 2;
t176 = 0.1e1 / t178 ^ 2;
t238 = t188 * t178;
t165 = -t259 * t179 + t187 * t238;
t160 = t165 ^ 2;
t141 = t144 * t160 + 0.1e1;
t180 = qJD(4) + qJD(5);
t215 = qJD(1) * t259;
t231 = qJD(3) * t189;
t201 = -t187 * t215 - t188 * t231;
t198 = t259 * t180 - t201;
t210 = t180 * t187 + qJD(1);
t149 = t198 * t178 + t210 * t237;
t253 = t149 * t144;
t164 = t167 ^ 2;
t184 = 0.1e1 / t189 ^ 2;
t243 = t176 * t184;
t156 = t164 * t243 + 0.1e1;
t153 = 0.1e1 / t156;
t233 = qJD(3) * t187;
t241 = t179 * t189;
t202 = -t178 * t233 + t180 * t241;
t222 = t167 * t243;
t218 = t259 * t189;
t207 = qJD(3) * t218;
t208 = t179 * t219;
t147 = -t180 * t208 - t178 * t207 - t179 * t215 + (qJD(1) * t187 + t180) * t238;
t244 = t175 * t183;
t225 = t147 * t244;
t135 = (-t202 * t222 - t225) * t153;
t200 = -t135 * t167 - t202;
t130 = (-t135 * t235 - t147) * t151 - t200 * t152;
t145 = t143 * t144;
t257 = t130 * t145;
t258 = (-t160 * t257 + t165 * t253) / t141 ^ 2;
t177 = t175 * t176;
t182 = t189 ^ 2;
t185 = t183 / t182;
t242 = t179 * t180;
t220 = t184 * t242;
t256 = (-t147 * t222 + (t176 * t185 * t233 - t177 * t220) * t164) / t156 ^ 2;
t181 = t188 ^ 2;
t240 = t181 * t182;
t224 = t162 * t240;
t159 = 0.1e1 + t224;
t199 = -t181 * t187 * t231 + t182 * t188 * t215;
t150 = t198 * t179 - t210 * t238;
t252 = t150 * t161 * t162;
t209 = t240 * t252;
t255 = (t199 * t162 - t209) / t159 ^ 2;
t254 = t144 * t165;
t251 = t151 * t165;
t249 = t151 * t189;
t248 = t152 * t165;
t247 = t152 * t167;
t246 = t152 * t187;
t245 = t167 * t180;
t239 = t187 * t188;
t234 = qJD(1) * t188;
t232 = qJD(3) * t188;
t230 = 0.2e1 * t258;
t229 = -0.2e1 * t256;
t228 = 0.2e1 * t255;
t227 = 0.2e1 * t145 * t165;
t226 = t144 * t251;
t223 = t167 * t244;
t221 = t175 * t184 * t187;
t217 = t184 * t233;
t204 = t167 * t221 + t259;
t142 = t204 * t153;
t216 = t259 - t142;
t214 = -0.2e1 * t143 * t258;
t213 = t144 * t230;
t212 = 0.2e1 * t183 * t256;
t211 = -0.2e1 * t165 * t236;
t206 = t175 * t212;
t168 = t208 - t238;
t205 = t167 * t176 * t179 - t168 * t175;
t203 = t162 * t168 * t188 - t259 * t161;
t197 = t149 * t244 - (t176 * t183 * t242 - t175 * t217) * t165;
t157 = 0.1e1 / t159;
t148 = t166 * qJD(1) - t179 * t207 + t245;
t139 = 0.1e1 / t141;
t138 = t205 * t183 * t153;
t134 = (-t151 + (-t152 * t223 + t151) * t153) * t165;
t133 = t142 * t247 + (t216 * t249 - t246) * t178;
t132 = t152 * t241 + t151 * t168 - (-t151 * t235 + t247) * t138;
t131 = t162 * t211 * t255 + (t211 * t252 + (t149 * t236 + (-t187 * t232 + t189 * t215) * t165) * t162) * t157;
t129 = t204 * t229 + (-t147 * t221 - t234 + (-t176 * t187 * t220 + (0.2e1 * t185 * t187 ^ 2 + t183) * t175 * qJD(3)) * t167) * t153;
t127 = t205 * t212 + (-t205 * t217 + ((-t148 + t245) * t175 + (0.2e1 * t167 * t177 * t242 + (-t168 * t180 + t147) * t176) * t179) * t183) * t153;
t126 = (t132 * t254 - t143 * t166) * t230 + (-t132 * t253 + t150 * t143 + (t132 * t227 - t144 * t166) * t130 - (-t179 * t233 - t180 * t235 + t127 * t167 + t138 * t147 + (t138 * t235 + t168) * t135) * t144 * t248 - (-t148 + (-t127 * t178 - t135 * t179) * t189 - t200 * t138) * t226) * t139;
t1 = [-t197 * t153 + t165 * t206, 0, t129, t127, t127, 0; t167 * t214 + (-t147 * t143 + (-t130 * t167 - t134 * t149) * t144) * t139 + (t134 * t213 + (0.2e1 * t134 * t257 + (-t149 * t153 + t149 - (t135 * t153 * t223 + t229) * t165) * t144 * t151 + (-(t167 * t206 - t135) * t254 + (-(t135 + t225) * t165 + t197 * t167) * t144 * t153) * t152) * t139) * t165, 0, t133 * t165 * t213 + (-(t129 * t247 + (-t135 * t250 - t147 * t152) * t142) * t254 + (t143 * t236 - (-t142 * t249 + t151 * t218 - t246) * t254) * t242 + (t130 * t227 - t253) * t133) * t139 + (t214 * t236 + ((-t143 * t232 - (-t216 * qJD(3) + t135) * t226) * t187 + (t143 * t215 + (-t188 * t130 - (-t129 - t234) * t251 - (t216 * t135 - qJD(3)) * t248) * t144) * t189) * t139) * t178, t126, t126, 0; t203 * t189 * t228 + (t203 * t233 + ((-qJD(1) * t161 + 0.2e1 * t168 * t252) * t188 + (t148 * t188 - t259 * t150 - t168 * t215) * t162) * t189) * t157, 0 (t161 * t239 + t179 * t224) * t228 + (0.2e1 * t179 * t209 + t201 * t161 + (t178 * t180 * t240 + t150 * t239 - 0.2e1 * t179 * t199) * t162) * t157, t131, t131, 0;];
JaD_rot  = t1;
