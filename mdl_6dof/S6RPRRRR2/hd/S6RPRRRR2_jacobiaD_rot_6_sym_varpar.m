% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:15:25
% EndTime: 2019-02-26 21:15:26
% DurationCPUTime: 0.86s
% Computational Cost: add. (5971->100), mult. (4025->203), div. (771->12), fcn. (4686->9), ass. (0->100)
t196 = qJ(3) + qJ(4);
t189 = sin(t196);
t182 = t189 ^ 2;
t191 = cos(t196);
t184 = 0.1e1 / t191 ^ 2;
t241 = t182 * t184;
t194 = qJ(1) + pkin(11);
t186 = sin(t194);
t260 = 0.2e1 * t186;
t259 = t189 * t241;
t187 = cos(t194);
t195 = qJ(5) + qJ(6);
t188 = sin(t195);
t190 = cos(t195);
t233 = t190 * t191;
t169 = t186 * t188 + t187 * t233;
t192 = qJD(5) + qJD(6);
t213 = t191 * t192 - qJD(1);
t193 = qJD(3) + qJD(4);
t234 = t189 * t193;
t258 = t213 * t188 + t190 * t234;
t238 = t186 * t189;
t173 = atan2(-t238, -t191);
t171 = cos(t173);
t170 = sin(t173);
t222 = t170 * t238;
t157 = -t171 * t191 - t222;
t154 = 0.1e1 / t157;
t163 = 0.1e1 / t169;
t183 = 0.1e1 / t191;
t155 = 0.1e1 / t157 ^ 2;
t164 = 0.1e1 / t169 ^ 2;
t257 = -0.2e1 * t189;
t179 = t186 ^ 2;
t176 = t179 * t241 + 0.1e1;
t174 = 0.1e1 / t176;
t256 = t174 - 0.1e1;
t228 = qJD(1) * t189;
t218 = t187 * t228;
t232 = t191 * t193;
t237 = t186 * t193;
t147 = (-(-t186 * t232 - t218) * t183 + t237 * t241) * t174;
t244 = t171 * t189;
t142 = (-t147 * t186 + t193) * t244 + (-t218 + (t147 - t237) * t191) * t170;
t255 = t142 * t154 * t155;
t254 = t147 * t170;
t253 = t147 * t189;
t212 = -qJD(1) * t191 + t192;
t208 = t190 * t212;
t149 = t186 * t208 - t258 * t187;
t252 = t149 * t163 * t164;
t205 = (t189 + t259) * t183 * t193;
t229 = qJD(1) * t187;
t210 = t182 * t186 * t229;
t251 = (t179 * t205 + t184 * t210) / t176 ^ 2;
t250 = t155 * t187;
t249 = t155 * t189;
t235 = t188 * t191;
t206 = t186 * t235 + t187 * t190;
t220 = t188 * t234;
t148 = t206 * qJD(1) - t169 * t192 + t187 * t220;
t168 = -t186 * t190 + t187 * t235;
t162 = t168 ^ 2;
t161 = t162 * t164 + 0.1e1;
t246 = t164 * t168;
t248 = 0.1e1 / t161 ^ 2 * (-t148 * t246 - t162 * t252);
t217 = 0.1e1 + t241;
t160 = t217 * t186 * t174;
t247 = t160 * t186;
t245 = t170 * t191;
t180 = t187 ^ 2;
t243 = t180 * t182;
t242 = t182 * t183;
t240 = t183 * t186;
t236 = t187 * t188;
t231 = t193 * t154;
t230 = qJD(1) * t186;
t152 = t155 * t243 + 0.1e1;
t227 = 0.2e1 * (-t243 * t255 + (t180 * t189 * t232 - t210) * t155) / t152 ^ 2;
t226 = 0.2e1 * t255;
t225 = 0.2e1 * t248;
t224 = t168 * t252;
t223 = t187 * t249;
t216 = t189 * t227;
t215 = t251 * t260;
t214 = t251 * t257;
t211 = t171 * t174 * t242;
t209 = t217 * t187;
t207 = -t163 * t188 + t190 * t246;
t204 = t207 * t189;
t167 = -t186 * t233 + t236;
t158 = 0.1e1 / t161;
t150 = 0.1e1 / t152;
t146 = (t256 * t189 * t170 - t186 * t211) * t187;
t145 = -t186 * t245 + t244 + (-t171 * t238 + t245) * t160;
t143 = -t217 * t215 + (qJD(1) * t209 + t205 * t260) * t174;
t140 = -0.2e1 * t248 + 0.2e1 * (-t148 * t164 * t158 + (-t158 * t252 - t164 * t248) * t168) * t168;
t139 = -t187 * t204 * t225 + (-t204 * t230 + (t207 * t232 + ((-t163 * t192 - 0.2e1 * t224) * t190 + (-t148 * t190 + (-t168 * t192 + t149) * t188) * t164) * t189) * t187) * t158;
t138 = (t145 * t249 - t154 * t191) * t187 * t227 + ((-t154 * t230 + (-t145 * t193 - t142) * t250) * t191 + (-t187 * t231 - (-t143 * t171 * t186 + t170 * t237 + t247 * t254 - t254 + (-t170 * t193 - t171 * t229) * t160) * t223 + (t155 * t230 + t187 * t226) * t145 - ((t143 - t229) * t170 + ((0.1e1 - t247) * t193 + (t160 - t186) * t147) * t171) * t191 * t250) * t189) * t150;
t1 = [t187 * t183 * t214 + (t193 * t209 - t228 * t240) * t174, 0, t143, t143, 0, 0; (t154 * t216 + (-t191 * t231 + (qJD(1) * t146 + t142) * t249) * t150) * t186 + (t155 * t216 * t146 + (-((t214 - t232 + (t147 * t182 * t240 + t232) * t174) * t170 + (t215 * t242 - t253 + (t253 + (t257 - t259) * t237) * t174) * t171) * t223 + (-t155 * t232 + t189 * t226) * t146 + (-t154 + ((-t179 + t180) * t211 + t256 * t222) * t155) * t228) * t150) * t187, 0, t138, t138, 0, 0; (t163 * t206 + t167 * t246) * t225 + (0.2e1 * t167 * t224 + (t167 * t148 + t206 * t149 + (-t258 * t186 - t187 * t208) * t168) * t164 + (t212 * t236 + (-t213 * t190 + t220) * t186) * t163) * t158, 0, t139, t139, t140, t140;];
JaD_rot  = t1;
