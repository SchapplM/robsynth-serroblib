% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:01
% EndTime: 2019-02-26 22:07:02
% DurationCPUTime: 0.99s
% Computational Cost: add. (1687->119), mult. (4276->253), div. (493->12), fcn. (5073->11), ass. (0->111)
t190 = sin(qJ(2));
t183 = t190 ^ 2;
t193 = cos(qJ(2));
t186 = 0.1e1 / t193 ^ 2;
t246 = t183 * t186;
t191 = sin(qJ(1));
t184 = t191 ^ 2;
t175 = t184 * t246 + 0.1e1;
t185 = 0.1e1 / t193;
t243 = t185 * t190;
t265 = t190 * t246;
t202 = qJD(2) * (t185 * t265 + t243);
t194 = cos(qJ(1));
t235 = qJD(1) * t194;
t244 = t183 * t191;
t211 = t235 * t244;
t249 = (t184 * t202 + t186 * t211) / t175 ^ 2;
t266 = -0.2e1 * t249;
t218 = 0.1e1 + t246;
t264 = t191 * t218;
t192 = cos(qJ(3));
t263 = (-qJD(3) + qJD(6)) * t192;
t189 = sin(qJ(3));
t231 = qJD(3) * t189;
t262 = -qJD(6) * t189 + t231;
t238 = t191 * t193;
t203 = t189 * t238 + t192 * t194;
t232 = qJD(2) * t194;
t219 = t190 * t232;
t237 = t193 * t194;
t221 = t192 * t237;
t140 = t203 * qJD(1) - qJD(3) * t221 + t189 * t219 - t191 * t231;
t213 = -qJD(1) * t193 + qJD(3);
t214 = qJD(3) * t193 - qJD(1);
t242 = t189 * t194;
t141 = -t214 * t242 + (t213 * t191 - t219) * t192;
t181 = pkin(10) + qJ(6);
t179 = sin(t181);
t180 = cos(t181);
t239 = t191 * t192;
t168 = t189 * t237 - t239;
t169 = t191 * t189 + t221;
t207 = t168 * t180 - t169 * t179;
t133 = t207 * qJD(6) - t140 * t179 + t141 * t180;
t153 = t168 * t179 + t169 * t180;
t145 = 0.1e1 / t153;
t205 = t179 * t189 + t180 * t192;
t206 = t179 * t192 - t180 * t189;
t146 = 0.1e1 / t153 ^ 2;
t252 = t146 * t207;
t261 = t206 * t145 + t205 * t252;
t240 = t191 * t190;
t174 = atan2(t240, t193);
t171 = cos(t174);
t170 = sin(t174);
t223 = t170 * t240;
t160 = t171 * t193 + t223;
t157 = 0.1e1 / t160;
t158 = 0.1e1 / t160 ^ 2;
t260 = 0.2e1 * t190;
t172 = 0.1e1 / t175;
t259 = t172 - 0.1e1;
t188 = t194 ^ 2;
t245 = t183 * t188;
t156 = t158 * t245 + 0.1e1;
t233 = qJD(2) * t193;
t220 = t190 * t235;
t234 = qJD(2) * t191;
t139 = ((t191 * t233 + t220) * t185 + t234 * t246) * t172;
t247 = t171 * t190;
t130 = (t139 * t191 - qJD(2)) * t247 + (t220 + (-t139 + t234) * t193) * t170;
t257 = t130 * t157 * t158;
t258 = (-t245 * t257 + (t188 * t190 * t233 - t211) * t158) / t156 ^ 2;
t147 = t145 * t146;
t256 = t133 * t147;
t132 = t153 * qJD(6) + t140 * t180 + t141 * t179;
t144 = t207 ^ 2;
t137 = t144 * t146 + 0.1e1;
t255 = 0.1e1 / t137 ^ 2 * (-t132 * t252 - t144 * t256);
t254 = t139 * t170;
t253 = t139 * t190;
t241 = t190 * t194;
t164 = t205 * t241;
t251 = t146 * t164;
t250 = t158 * t190;
t162 = t172 * t264;
t248 = t162 * t191;
t236 = qJD(1) * t191;
t227 = -0.2e1 * t257;
t226 = 0.2e1 * t255;
t225 = -0.2e1 * t147 * t207;
t224 = t158 * t241;
t222 = t172 * t183 * t185;
t217 = -0.2e1 * t190 * t258;
t216 = t133 * t225;
t215 = t185 * t266;
t212 = t191 * t222;
t210 = t218 * t194;
t167 = -t192 * t238 + t242;
t208 = -t167 * t179 - t180 * t203;
t149 = t167 * t180 - t179 * t203;
t204 = t213 * t194;
t163 = t206 * t241;
t154 = 0.1e1 / t156;
t143 = t192 * t204 + (qJD(2) * t190 * t192 + t214 * t189) * t191;
t142 = -t214 * t239 + (t190 * t234 + t204) * t189;
t138 = (-t259 * t190 * t170 + t171 * t212) * t194;
t135 = 0.1e1 / t137;
t134 = t170 * t238 - t247 + (-t170 * t193 + t171 * t240) * t162;
t131 = t264 * t266 + (qJD(1) * t210 + 0.2e1 * t191 * t202) * t172;
t1 = [t215 * t241 + (qJD(2) * t210 - t236 * t243) * t172, t131, 0, 0, 0, 0; (t157 * t217 + (t157 * t233 + (-qJD(1) * t138 - t130) * t250) * t154) * t191 + (t158 * t217 * t138 + (((-t139 * t212 - t259 * t233 + t249 * t260) * t170 + (t215 * t244 + t253 + (-t253 + (t260 + t265) * t234) * t172) * t171) * t224 + (t158 * t233 + t190 * t227) * t138 + (t157 + ((-t184 + t188) * t171 * t222 + t259 * t223) * t158) * t190 * qJD(1)) * t154) * t194, 0.2e1 * (-t134 * t250 + t157 * t193) * t194 * t258 + ((t157 * t236 + (qJD(2) * t134 + t130) * t194 * t158) * t193 + (t157 * t232 + (t131 * t171 * t191 - t170 * t234 - t248 * t254 + t254 + (qJD(2) * t170 + t171 * t235) * t162) * t224 + (-t158 * t236 + t194 * t227) * t134 + ((-t131 + t235) * t170 + ((-0.1e1 + t248) * qJD(2) + (-t162 + t191) * t139) * t171) * t158 * t237) * t190) * t154, 0, 0, 0, 0; (t145 * t208 - t149 * t252) * t226 + ((t149 * qJD(6) - t142 * t180 + t143 * t179) * t145 + t149 * t216 + (t208 * t133 + (t208 * qJD(6) + t142 * t179 + t143 * t180) * t207 - t149 * t132) * t146) * t135 (t145 * t163 + t207 * t251) * t226 + (t132 * t251 + (t146 * t163 - t164 * t225) * t133 - t261 * t193 * t232 + (t261 * t236 + ((-t263 * t145 + t262 * t252) * t180 + (t262 * t145 + t263 * t252) * t179) * t194) * t190) * t135 (t145 * t153 + t207 * t252) * t226 + (-t133 * t145 - t207 * t216 + (0.2e1 * t207 * t132 + t153 * t133) * t146) * t135, 0, 0, -0.2e1 * t255 - 0.2e1 * (t132 * t146 * t135 - (-t135 * t256 - t146 * t255) * t207) * t207;];
JaD_rot  = t1;
