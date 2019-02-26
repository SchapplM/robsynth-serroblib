% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRP5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:11:36
% EndTime: 2019-02-26 22:11:37
% DurationCPUTime: 1.08s
% Computational Cost: add. (11090->123), mult. (8378->264), div. (1558->15), fcn. (10537->9), ass. (0->116)
t187 = qJ(3) + pkin(10) + qJ(5);
t186 = cos(t187);
t270 = 0.2e1 * t186;
t185 = sin(t187);
t195 = cos(qJ(2));
t196 = cos(qJ(1));
t243 = t196 * t186;
t226 = t195 * t243;
t264 = sin(qJ(1));
t171 = t185 * t264 + t226;
t165 = 0.1e1 / t171 ^ 2;
t194 = sin(qJ(2));
t189 = t194 ^ 2;
t193 = t196 ^ 2;
t248 = t189 * t193;
t227 = t165 * t248;
t161 = 0.1e1 + t227;
t219 = qJD(1) * t264;
t241 = qJD(2) * t195;
t204 = t189 * t196 * t219 - t193 * t194 * t241;
t188 = qJD(3) + qJD(5);
t240 = qJD(2) * t196;
t221 = t194 * t240;
t207 = t195 * t219 + t221;
t225 = t264 * t188;
t244 = t196 * t185;
t150 = (-t188 * t195 + qJD(1)) * t244 + (t225 - t207) * t186;
t164 = 0.1e1 / t171;
t259 = t150 * t164 * t165;
t214 = t248 * t259;
t269 = (-t165 * t204 - t214) / t161 ^ 2;
t223 = t264 * t195;
t167 = t185 * t223 + t243;
t268 = t167 * t188;
t246 = t194 * t196;
t149 = qJD(1) * t167 - t188 * t226 + (t221 - t225) * t185;
t170 = -t186 * t264 + t195 * t244;
t182 = 0.1e1 / t185;
t183 = 0.1e1 / t185 ^ 2;
t190 = 0.1e1 / t194;
t191 = 0.1e1 / t194 ^ 2;
t222 = t191 * t241;
t250 = t186 * t188;
t252 = t182 * t190;
t267 = t170 * (t183 * t190 * t250 + t182 * t222) + t149 * t252;
t247 = t194 * t185;
t157 = atan2(-t167, t247);
t154 = cos(t157);
t153 = sin(t157);
t258 = t153 * t167;
t148 = t154 * t247 - t258;
t145 = 0.1e1 / t148;
t146 = 0.1e1 / t148 ^ 2;
t266 = -0.2e1 * t167;
t265 = 0.2e1 * t170;
t162 = t167 ^ 2;
t251 = t183 * t191;
t158 = t162 * t251 + 0.1e1;
t155 = 0.1e1 / t158;
t249 = t186 * t194;
t208 = t185 * t241 + t188 * t249;
t230 = t167 * t251;
t224 = t264 * t194;
t212 = qJD(2) * t224;
t242 = qJD(1) * t196;
t151 = -t185 * t212 - t188 * t244 - t186 * t219 + (t185 * t242 + t186 * t225) * t195;
t232 = t151 * t252;
t137 = (t208 * t230 - t232) * t155;
t205 = -t137 * t167 + t208;
t132 = (-t137 * t247 - t151) * t153 + t205 * t154;
t147 = t145 * t146;
t263 = t132 * t147;
t184 = t182 * t183;
t192 = t190 / t189;
t228 = t191 * t250;
t262 = (t151 * t230 + (-t183 * t192 * t241 - t184 * t228) * t162) / t158 ^ 2;
t261 = t146 * t170;
t260 = t149 * t146;
t257 = t153 * t170;
t256 = t153 * t194;
t255 = t154 * t167;
t254 = t154 * t170;
t253 = t154 * t195;
t245 = t195 * t196;
t163 = t170 ^ 2;
t143 = t163 * t146 + 0.1e1;
t239 = 0.2e1 * (-t163 * t263 - t170 * t260) / t143 ^ 2;
t238 = -0.2e1 * t262;
t237 = 0.2e1 * t269;
t236 = t147 * t265;
t235 = t190 * t262;
t234 = t146 * t257;
t231 = t167 * t252;
t229 = t182 * t191 * t195;
t210 = t167 * t229 + t264;
t144 = t210 * t155;
t220 = t264 - t144;
t218 = t145 * t239;
t217 = t146 * t239;
t216 = t246 * t265;
t215 = t182 * t235;
t169 = t186 * t223 - t244;
t211 = t167 * t183 * t186 - t169 * t182;
t209 = t165 * t169 * t196 - t164 * t264;
t159 = 0.1e1 / t161;
t152 = qJD(1) * t171 - t186 * t212 - t268;
t141 = 0.1e1 / t143;
t140 = t211 * t190 * t155;
t136 = (-t153 + (t154 * t231 + t153) * t155) * t170;
t135 = -t144 * t255 + (t220 * t256 + t253) * t185;
t134 = t154 * t249 - t153 * t169 + (-t153 * t247 - t255) * t140;
t133 = t165 * t216 * t269 + (t216 * t259 + (t149 * t246 + (t194 * t219 - t195 * t240) * t170) * t165) * t159;
t131 = t210 * t238 + (t151 * t229 + t242 + (-t183 * t195 * t228 + (-0.2e1 * t192 * t195 ^ 2 - t190) * t182 * qJD(2)) * t167) * t155;
t129 = -0.2e1 * t211 * t235 + (-t211 * t222 + ((-t152 - t268) * t182 + (t184 * t250 * t266 + (t169 * t188 + t151) * t183) * t186) * t190) * t155;
t128 = (t134 * t261 - t145 * t171) * t239 + (t134 * t260 + t150 * t145 + (t134 * t236 - t171 * t146) * t132 - (t186 * t241 - t188 * t247 - t129 * t167 - t140 * t151 + (-t140 * t247 - t169) * t137) * t146 * t254 - (-t152 + (-t129 * t185 - t137 * t186) * t194 - t205 * t140) * t234) * t141;
t1 = [t267 * t155 + t215 * t265, t131, t129, 0, t129, 0; t167 * t218 + (-t151 * t145 + (t132 * t167 + t136 * t149) * t146) * t141 + (t136 * t217 + (0.2e1 * t136 * t263 + (t149 * t155 - t149 - (-t137 * t155 * t231 + t238) * t170) * t146 * t153 + (-(t215 * t266 - t137) * t261 + (-(t137 + t232) * t170 + t267 * t167) * t146 * t155) * t154) * t141) * t170, t135 * t170 * t217 + (-(-t131 * t255 + (t137 * t258 - t151 * t154) * t144) * t261 + (-t145 * t246 - (-t144 * t256 + t153 * t224 + t253) * t261) * t250 + (t236 * t132 + t260) * t135) * t141 + (t218 * t246 + ((-t145 * t240 - (qJD(2) * t220 - t137) * t234) * t195 + (t145 * t219 + (t196 * t132 - (-t131 + t242) * t257 - (t137 * t220 - qJD(2)) * t254) * t146) * t194) * t141) * t185, t128, 0, t128, 0; t209 * t194 * t237 + (-t209 * t241 + ((qJD(1) * t164 + 0.2e1 * t169 * t259) * t196 + (-t150 * t264 - t152 * t196 + t169 * t219) * t165) * t194) * t159 (t164 * t245 + t186 * t227) * t237 + (t214 * t270 + t207 * t164 + (t185 * t188 * t248 + t150 * t245 + t204 * t270) * t165) * t159, t133, 0, t133, 0;];
JaD_rot  = t1;
