% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:09:03
% EndTime: 2019-02-26 21:09:04
% DurationCPUTime: 1.09s
% Computational Cost: add. (9974->127), mult. (8378->272), div. (1558->15), fcn. (10537->9), ass. (0->117)
t192 = qJ(4) + qJ(5);
t185 = sin(t192);
t186 = cos(t192);
t240 = qJ(1) + pkin(10);
t222 = sin(t240);
t184 = cos(t240);
t194 = cos(qJ(3));
t249 = t184 * t194;
t169 = t222 * t185 + t186 * t249;
t163 = 0.1e1 / t169 ^ 2;
t180 = t184 ^ 2;
t193 = sin(qJ(3));
t188 = t193 ^ 2;
t253 = t180 * t188;
t231 = t163 * t253;
t158 = 0.1e1 + t231;
t212 = qJD(1) * t222;
t208 = t194 * t212;
t187 = qJD(4) + qJD(5);
t216 = t222 * t187;
t242 = qJD(3) * t193;
t245 = t187 * t194;
t148 = (-t208 + t216) * t186 + (-t186 * t242 + (qJD(1) - t245) * t185) * t184;
t162 = 0.1e1 / t169;
t260 = t148 * t162 * t163;
t217 = t253 * t260;
t241 = qJD(3) * t194;
t223 = t193 * t241;
t269 = (-t217 + (-t184 * t188 * t212 + t180 * t223) * t163) / t158 ^ 2;
t182 = 0.1e1 / t185 ^ 2;
t247 = t186 * t187;
t268 = t182 * t247;
t250 = t184 * t193;
t214 = t222 * t194;
t211 = t185 * t214;
t165 = t184 * t186 + t211;
t225 = t184 * t242;
t226 = t186 * t245;
t147 = t165 * qJD(1) - t184 * t226 + (-t216 + t225) * t185;
t168 = t185 * t249 - t222 * t186;
t181 = 0.1e1 / t185;
t189 = 0.1e1 / t193;
t190 = 0.1e1 / t193 ^ 2;
t224 = t190 * t241;
t252 = t181 * t189;
t267 = (t181 * t224 + t189 * t268) * t168 + t147 * t252;
t244 = t193 * t185;
t157 = atan2(-t165, t244);
t152 = cos(t157);
t151 = sin(t157);
t259 = t151 * t165;
t146 = t152 * t244 - t259;
t143 = 0.1e1 / t146;
t144 = 0.1e1 / t146 ^ 2;
t266 = -0.2e1 * t165;
t265 = 0.2e1 * t168;
t160 = t165 ^ 2;
t251 = t182 * t190;
t159 = t160 * t251 + 0.1e1;
t155 = 0.1e1 / t159;
t246 = t186 * t193;
t206 = t185 * t241 + t187 * t246;
t229 = t165 * t251;
t215 = t222 * t193;
t209 = qJD(3) * t215;
t210 = t186 * t212;
t243 = qJD(1) * t184;
t248 = t185 * t187;
t149 = -t185 * t209 - t184 * t248 - t210 + (t185 * t243 + t186 * t216) * t194;
t232 = t149 * t252;
t135 = (t206 * t229 - t232) * t155;
t202 = -t135 * t165 + t206;
t131 = (-t135 * t244 - t149) * t151 + t202 * t152;
t145 = t143 * t144;
t264 = t131 * t145;
t191 = t189 / t188;
t227 = t181 * t268;
t263 = (t149 * t229 + (-t182 * t191 * t241 - t190 * t227) * t160) / t159 ^ 2;
t262 = t144 * t168;
t261 = t147 * t144;
t258 = t151 * t168;
t257 = t151 * t193;
t256 = t152 * t165;
t255 = t152 * t168;
t254 = t152 * t194;
t161 = t168 ^ 2;
t141 = t144 * t161 + 0.1e1;
t239 = 0.2e1 * (-t161 * t264 - t168 * t261) / t141 ^ 2;
t238 = -0.2e1 * t263;
t237 = 0.2e1 * t269;
t236 = t145 * t265;
t235 = t189 * t263;
t234 = t144 * t258;
t230 = t165 * t252;
t228 = t181 * t190 * t194;
t221 = t143 * t239;
t220 = t144 * t239;
t219 = t250 * t265;
t218 = t181 * t235;
t205 = t165 * t228 + t222;
t142 = t205 * t155;
t213 = t222 - t142;
t167 = -t184 * t185 + t186 * t214;
t207 = t165 * t182 * t186 - t167 * t181;
t204 = t163 * t167 * t184 - t222 * t162;
t153 = 0.1e1 / t158;
t150 = t169 * qJD(1) - t184 * t247 - t186 * t209 - t187 * t211;
t139 = 0.1e1 / t141;
t138 = t207 * t189 * t155;
t134 = (-t151 + (t152 * t230 + t151) * t155) * t168;
t133 = -t142 * t256 + (t213 * t257 + t254) * t185;
t132 = t152 * t246 - t151 * t167 + (-t151 * t244 - t256) * t138;
t130 = t163 * t219 * t269 + (t219 * t260 + (t147 * t250 + (-t184 * t241 + t193 * t212) * t168) * t163) * t153;
t129 = t205 * t238 + (t149 * t228 + t243 + (-t226 * t251 + (-0.2e1 * t191 * t194 ^ 2 - t189) * t181 * qJD(3)) * t165) * t155;
t127 = -0.2e1 * t207 * t235 + (-t207 * t224 + ((-t165 * t187 - t150) * t181 + (t227 * t266 + (t167 * t187 + t149) * t182) * t186) * t189) * t155;
t126 = (t132 * t262 - t143 * t169) * t239 + (t132 * t261 + t148 * t143 + (t132 * t236 - t144 * t169) * t131 - (t186 * t241 - t187 * t244 - t127 * t165 - t138 * t149 + (-t138 * t244 - t167) * t135) * t144 * t255 - (-t150 + (-t127 * t185 - t135 * t186) * t193 - t202 * t138) * t234) * t139;
t1 = [t267 * t155 + t218 * t265, 0, t129, t127, t127, 0; t165 * t221 + (-t149 * t143 + (t131 * t165 + t134 * t147) * t144) * t139 + (t134 * t220 + (0.2e1 * t134 * t264 + (t147 * t155 - t147 - (-t135 * t155 * t230 + t238) * t168) * t144 * t151 + (-(t218 * t266 - t135) * t262 + (-(t135 + t232) * t168 + t267 * t165) * t144 * t155) * t152) * t139) * t168, 0, t133 * t168 * t220 + (-(-t129 * t256 + (t135 * t259 - t149 * t152) * t142) * t262 + (-t143 * t250 - (-t142 * t257 + t151 * t215 + t254) * t262) * t247 + (t131 * t236 + t261) * t133) * t139 + (t221 * t250 + ((-t184 * qJD(3) * t143 - (t213 * qJD(3) - t135) * t234) * t194 + (t143 * t212 + (t184 * t131 - (-t129 + t243) * t258 - (t213 * t135 - qJD(3)) * t255) * t144) * t193) * t139) * t185, t126, t126, 0; t204 * t193 * t237 + (-t204 * t241 + ((qJD(1) * t162 + 0.2e1 * t167 * t260) * t184 + (-t222 * t148 - t150 * t184 + t167 * t212) * t163) * t193) * t153, 0 (t162 * t249 + t186 * t231) * t237 + (0.2e1 * t186 * t217 + (t208 + t225) * t162 + ((t148 * t194 + 0.2e1 * t188 * t210) * t184 + (-0.2e1 * t186 * t223 + t188 * t248) * t180) * t163) * t153, t130, t130, 0;];
JaD_rot  = t1;
