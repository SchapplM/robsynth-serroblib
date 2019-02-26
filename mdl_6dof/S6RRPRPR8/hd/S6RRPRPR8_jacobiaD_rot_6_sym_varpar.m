% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:41:47
% EndTime: 2019-02-26 21:41:48
% DurationCPUTime: 1.07s
% Computational Cost: add. (2275->116), mult. (4276->247), div. (493->12), fcn. (5073->11), ass. (0->112)
t269 = qJD(4) - qJD(6);
t190 = sin(qJ(2));
t183 = t190 ^ 2;
t193 = cos(qJ(2));
t186 = 0.1e1 / t193 ^ 2;
t247 = t183 * t186;
t191 = sin(qJ(1));
t184 = t191 ^ 2;
t176 = t184 * t247 + 0.1e1;
t185 = 0.1e1 / t193;
t244 = t185 * t190;
t267 = t190 * t247;
t203 = qJD(2) * (t185 * t267 + t244);
t194 = cos(qJ(1));
t235 = qJD(1) * t194;
t245 = t183 * t191;
t211 = t235 * t245;
t250 = (t184 * t203 + t186 * t211) / t176 ^ 2;
t268 = -0.2e1 * t250;
t218 = 0.1e1 + t247;
t266 = t191 * t218;
t181 = pkin(10) + qJ(4);
t180 = cos(t181);
t237 = t194 * t180;
t179 = sin(t181);
t242 = t191 * t179;
t169 = t193 * t237 + t242;
t192 = cos(qJ(6));
t265 = t269 * t192;
t189 = sin(qJ(6));
t264 = t269 * t189;
t239 = t191 * t193;
t204 = t179 * t239 + t237;
t232 = qJD(2) * t194;
t219 = t190 * t232;
t139 = t204 * qJD(1) - t169 * qJD(4) + t179 * t219;
t213 = -qJD(1) * t193 + qJD(4);
t214 = qJD(4) * t193 - qJD(1);
t238 = t194 * t179;
t140 = -t214 * t238 + (t213 * t191 - t219) * t180;
t241 = t191 * t180;
t168 = t193 * t238 - t241;
t207 = t168 * t192 - t169 * t189;
t133 = t207 * qJD(6) - t139 * t189 + t140 * t192;
t153 = t168 * t189 + t169 * t192;
t145 = 0.1e1 / t153;
t205 = t179 * t189 + t180 * t192;
t206 = t179 * t192 - t180 * t189;
t146 = 0.1e1 / t153 ^ 2;
t254 = t146 * t207;
t263 = t206 * t145 - t205 * t254;
t240 = t191 * t190;
t175 = atan2(t240, t193);
t172 = cos(t175);
t171 = sin(t175);
t223 = t171 * t240;
t160 = t172 * t193 + t223;
t157 = 0.1e1 / t160;
t158 = 0.1e1 / t160 ^ 2;
t262 = 0.2e1 * t190;
t173 = 0.1e1 / t176;
t261 = t173 - 0.1e1;
t188 = t194 ^ 2;
t246 = t183 * t188;
t156 = t158 * t246 + 0.1e1;
t233 = qJD(2) * t193;
t220 = t190 * t235;
t234 = qJD(2) * t191;
t143 = ((t191 * t233 + t220) * t185 + t234 * t247) * t173;
t248 = t172 * t190;
t130 = (t143 * t191 - qJD(2)) * t248 + (t220 + (-t143 + t234) * t193) * t171;
t259 = t130 * t157 * t158;
t260 = (-t246 * t259 + (t188 * t190 * t233 - t211) * t158) / t156 ^ 2;
t147 = t145 * t146;
t258 = t133 * t147;
t132 = t153 * qJD(6) + t139 * t192 + t140 * t189;
t144 = t207 ^ 2;
t137 = t144 * t146 + 0.1e1;
t257 = 0.1e1 / t137 ^ 2 * (-t132 * t254 - t144 * t258);
t256 = t143 * t171;
t255 = t143 * t190;
t243 = t190 * t194;
t164 = t205 * t243;
t253 = t146 * t164;
t252 = t158 * t190;
t251 = t158 * t194;
t162 = t173 * t266;
t249 = t162 * t191;
t236 = qJD(1) * t191;
t227 = -0.2e1 * t259;
t226 = 0.2e1 * t257;
t225 = -0.2e1 * t147 * t207;
t224 = t158 * t243;
t222 = t173 * t183 * t185;
t217 = -0.2e1 * t190 * t260;
t216 = t133 * t225;
t215 = t185 * t268;
t212 = t191 * t222;
t210 = t194 * t218;
t167 = -t180 * t239 + t238;
t208 = -t167 * t189 - t192 * t204;
t149 = t167 * t192 - t189 * t204;
t202 = t190 * t234 + t213 * t194;
t163 = t206 * t243;
t154 = 0.1e1 / t156;
t142 = t202 * t180 + t214 * t242;
t141 = t202 * t179 - t214 * t241;
t138 = (-t261 * t190 * t171 + t172 * t212) * t194;
t135 = 0.1e1 / t137;
t134 = t171 * t239 - t248 + (-t171 * t193 + t172 * t240) * t162;
t131 = t266 * t268 + (qJD(1) * t210 + 0.2e1 * t191 * t203) * t173;
t1 = [t215 * t243 + (qJD(2) * t210 - t236 * t244) * t173, t131, 0, 0, 0, 0; (t157 * t217 + (t157 * t233 + (-qJD(1) * t138 - t130) * t252) * t154) * t191 + (t158 * t217 * t138 + (((-t143 * t212 - t261 * t233 + t250 * t262) * t171 + (t215 * t245 + t255 + (-t255 + (t262 + t267) * t234) * t173) * t172) * t224 + (t158 * t233 + t190 * t227) * t138 + (t157 + ((-t184 + t188) * t172 * t222 + t261 * t223) * t158) * t190 * qJD(1)) * t154) * t194, 0.2e1 * (-t134 * t252 + t157 * t193) * t194 * t260 + ((t157 * t236 + (qJD(2) * t134 + t130) * t251) * t193 + (t157 * t232 + (t131 * t172 * t191 - t171 * t234 - t249 * t256 + t256 + (qJD(2) * t171 + t172 * t235) * t162) * t224 + (-t158 * t236 + t194 * t227) * t134 + ((-t131 + t235) * t171 + ((-0.1e1 + t249) * qJD(2) + (-t162 + t191) * t143) * t172) * t193 * t251) * t190) * t154, 0, 0, 0, 0; (t145 * t208 - t149 * t254) * t226 + ((t149 * qJD(6) - t141 * t192 + t142 * t189) * t145 + t149 * t216 + (t208 * t133 + (t208 * qJD(6) + t141 * t189 + t142 * t192) * t207 - t149 * t132) * t146) * t135 (-t145 * t163 + t207 * t253) * t226 + (t132 * t253 + (-t163 * t146 - t164 * t225) * t133 + t263 * t193 * t232 + (-t263 * t236 + ((t265 * t145 - t264 * t254) * t180 + (t264 * t145 + t265 * t254) * t179) * t194) * t190) * t135, 0 (t145 * t153 + t207 * t254) * t226 + (-t133 * t145 - t207 * t216 + (0.2e1 * t207 * t132 + t153 * t133) * t146) * t135, 0, -0.2e1 * t257 - 0.2e1 * (t132 * t146 * t135 - (-t135 * t258 - t146 * t257) * t207) * t207;];
JaD_rot  = t1;
