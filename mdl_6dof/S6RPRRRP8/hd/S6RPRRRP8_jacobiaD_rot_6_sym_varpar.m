% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP8
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
% Datum: 2019-02-26 21:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:12:04
% EndTime: 2019-02-26 21:12:05
% DurationCPUTime: 1.04s
% Computational Cost: add. (5529->123), mult. (8382->269), div. (1515->15), fcn. (10508->9), ass. (0->119)
t188 = qJ(3) + qJ(4);
t182 = cos(t188);
t190 = sin(qJ(1));
t246 = t182 * t190;
t181 = sin(t188);
t189 = sin(qJ(5));
t265 = cos(qJ(1));
t221 = t265 * t189;
t191 = cos(qJ(5));
t239 = t190 * t191;
t169 = t181 * t221 + t239;
t247 = t182 * t189;
t159 = atan2(t169, t247);
t154 = cos(t159);
t153 = sin(t159);
t256 = t153 * t169;
t148 = t154 * t247 + t256;
t145 = 0.1e1 / t148;
t168 = t181 * t239 + t221;
t163 = 0.1e1 / t168;
t178 = 0.1e1 / t182;
t184 = 0.1e1 / t189;
t146 = 0.1e1 / t148 ^ 2;
t164 = 0.1e1 / t168 ^ 2;
t185 = 0.1e1 / t189 ^ 2;
t220 = t265 * t191;
t240 = t190 * t189;
t167 = t181 * t240 - t220;
t162 = t167 ^ 2;
t143 = t162 * t146 + 0.1e1;
t217 = qJD(1) * t265;
t207 = t181 * t217;
t201 = t265 * qJD(5) + t207;
t212 = qJD(5) * t181 + qJD(1);
t183 = qJD(3) + qJD(4);
t244 = t183 * t190;
t225 = t182 * t244;
t151 = t212 * t239 + (t201 + t225) * t189;
t259 = t151 * t146;
t166 = t169 ^ 2;
t179 = 0.1e1 / t182 ^ 2;
t249 = t179 * t185;
t161 = t166 * t249 + 0.1e1;
t157 = 0.1e1 / t161;
t236 = qJD(5) * t191;
t248 = t181 * t183;
t202 = t182 * t236 - t189 * t248;
t228 = t169 * t249;
t208 = t191 * t217;
t209 = t181 * t220;
t222 = t182 * t265;
t210 = t183 * t222;
t149 = -t189 * t210 - qJD(5) * t209 - t208 + (qJD(1) * t181 + qJD(5)) * t240;
t250 = t178 * t184;
t230 = t149 * t250;
t137 = (-t202 * t228 - t230) * t157;
t200 = t137 * t169 + t202;
t133 = (-t137 * t247 - t149) * t153 + t200 * t154;
t147 = t145 * t146;
t263 = t133 * t147;
t264 = (-t162 * t263 + t167 * t259) / t143 ^ 2;
t177 = t182 ^ 2;
t187 = t190 ^ 2;
t251 = t177 * t187;
t223 = t164 * t251;
t160 = 0.1e1 + t223;
t243 = t183 * t191;
t224 = t182 * t243;
t152 = t201 * t191 + (-t212 * t189 + t224) * t190;
t258 = t152 * t163 * t164;
t211 = t251 * t258;
t262 = (-t211 + (t177 * t190 * t217 - t182 * t187 * t248) * t164) / t160 ^ 2;
t180 = t178 / t177;
t186 = t184 * t185;
t261 = (-t149 * t228 + (-t179 * t186 * t236 + t180 * t185 * t248) * t166) / t161 ^ 2;
t260 = t146 * t167;
t257 = t153 * t167;
t255 = t153 * t182;
t254 = t154 * t167;
t253 = t154 * t169;
t252 = t154 * t181;
t245 = t183 * t184;
t242 = t185 * t191;
t241 = t190 * t145;
t238 = qJD(1) * t190;
t237 = qJD(5) * t189;
t235 = 0.2e1 * t264;
t234 = 0.2e1 * t262;
t233 = -0.2e1 * t261;
t232 = 0.2e1 * t147 * t167;
t231 = t146 * t257;
t229 = t169 * t250;
t227 = t179 * t181 * t184;
t226 = t181 * t245;
t219 = t185 * t236;
t204 = t169 * t227 + t265;
t144 = t204 * t157;
t218 = t265 - t144;
t216 = -0.2e1 * t145 * t264;
t215 = t146 * t235;
t214 = 0.2e1 * t178 * t261;
t213 = -0.2e1 * t167 * t246;
t206 = t184 * t214;
t170 = t209 - t240;
t205 = t169 * t242 - t170 * t184;
t203 = t164 * t170 * t190 - t265 * t163;
t199 = t151 * t250 - (t178 * t219 - t179 * t226) * t167;
t155 = 0.1e1 / t160;
t150 = t168 * qJD(1) + t169 * qJD(5) - t191 * t210;
t141 = 0.1e1 / t143;
t140 = t205 * t178 * t157;
t136 = (-t153 + (-t154 * t229 + t153) * t157) * t167;
t135 = t144 * t253 + (t218 * t255 - t252) * t189;
t134 = t154 * t182 * t191 + t153 * t170 - (-t153 * t247 + t253) * t140;
t132 = t204 * t233 + (-t149 * t227 - t238 + (t178 * t245 + (-t179 * t219 + 0.2e1 * t180 * t226) * t181) * t169) * t157;
t130 = t205 * t214 + (-t205 * t179 * t248 + (t149 * t242 - t150 * t184 + (-t170 * t242 + (0.2e1 * t186 * t191 ^ 2 + t184) * t169) * qJD(5)) * t178) * t157;
t129 = (t163 * t181 * t190 + t191 * t223) * t234 + (0.2e1 * t191 * t211 + (-t207 - t225) * t163 + ((t152 * t190 + 0.2e1 * t187 * t224) * t181 + (t187 * t237 - 0.2e1 * t190 * t208) * t177) * t164) * t155;
t128 = t135 * t167 * t215 + (-(t132 * t253 + (-t137 * t256 - t149 * t154) * t144) * t260 + (t133 * t232 - t259) * t135 + (t182 * t241 - (-t144 * t255 + t153 * t222 - t252) * t260) * t236) * t141 + (t216 * t246 + ((-t183 * t241 - (-t218 * t183 + t137) * t231) * t181 + (t145 * t217 + (-t190 * t133 - (-t132 - t238) * t257 - (t218 * t137 - t183) * t254) * t146) * t182) * t141) * t189;
t1 = [-t199 * t157 + t167 * t206, 0, t132, t132, t130, 0; t169 * t216 + (-t149 * t145 + (-t133 * t169 - t136 * t151) * t146) * t141 + (t136 * t215 + (0.2e1 * t136 * t263 + (-t151 * t157 + t151 - (t137 * t157 * t229 + t233) * t167) * t146 * t153 + (-(t169 * t206 - t137) * t260 + (-(t137 + t230) * t167 + t199 * t169) * t146 * t157) * t154) * t141) * t167, 0, t128, t128 (t134 * t260 - t145 * t168) * t235 + (-t134 * t259 + t152 * t145 + (t134 * t232 - t168 * t146) * t133 - (-t182 * t237 - t181 * t243 + t130 * t169 + t140 * t149 + (t140 * t247 + t170) * t137) * t146 * t254 - (-t150 + (-t130 * t189 - t137 * t191) * t182 + t200 * t140) * t231) * t141, 0; t203 * t182 * t234 + (t203 * t248 + ((-qJD(1) * t163 + 0.2e1 * t170 * t258) * t190 + (t150 * t190 - t265 * t152 - t170 * t217) * t164) * t182) * t155, 0, t129, t129, t164 * t213 * t262 + (t213 * t258 + (t151 * t246 + (-t181 * t244 + t182 * t217) * t167) * t164) * t155, 0;];
JaD_rot  = t1;
