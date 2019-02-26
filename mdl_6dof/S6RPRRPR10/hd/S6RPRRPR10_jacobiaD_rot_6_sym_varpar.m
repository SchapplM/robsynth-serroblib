% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR10_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:04
% EndTime: 2019-02-26 21:06:05
% DurationCPUTime: 1.06s
% Computational Cost: add. (1454->116), mult. (4276->251), div. (493->12), fcn. (5073->11), ass. (0->111)
t264 = -qJD(4) + qJD(6);
t186 = sin(qJ(3));
t178 = 0.1e1 / t186 ^ 2;
t190 = cos(qJ(3));
t182 = t190 ^ 2;
t244 = t178 * t182;
t263 = t190 * t244;
t191 = cos(qJ(1));
t217 = 0.1e1 + t244;
t262 = t191 * t217;
t189 = cos(qJ(4));
t261 = t264 * t189;
t185 = sin(qJ(4));
t260 = t264 * t185;
t187 = sin(qJ(1));
t213 = qJD(1) * t186 + qJD(4);
t201 = t213 * t191;
t214 = qJD(4) * t186 + qJD(1);
t202 = t214 * t189;
t232 = qJD(3) * t190;
t142 = t187 * t202 + (t187 * t232 + t201) * t185;
t203 = t214 * t185;
t143 = t189 * t201 + (t189 * t232 - t203) * t187;
t184 = sin(qJ(6));
t188 = cos(qJ(6));
t238 = t191 * t189;
t241 = t187 * t185;
t166 = t186 * t241 - t238;
t240 = t187 * t189;
t167 = t185 * t191 + t186 * t240;
t207 = t166 * t188 - t167 * t184;
t133 = qJD(6) * t207 + t142 * t184 + t143 * t188;
t154 = t166 * t184 + t167 * t188;
t148 = 0.1e1 / t154;
t204 = t184 * t185 + t188 * t189;
t205 = t184 * t189 - t185 * t188;
t149 = 0.1e1 / t154 ^ 2;
t251 = t149 * t207;
t259 = t148 * t205 + t204 * t251;
t237 = t191 * t190;
t172 = atan2(t237, -t186);
t170 = sin(t172);
t171 = cos(t172);
t160 = t170 * t237 - t171 * t186;
t157 = 0.1e1 / t160;
t177 = 0.1e1 / t186;
t158 = 0.1e1 / t160 ^ 2;
t183 = t191 ^ 2;
t175 = t183 * t244 + 0.1e1;
t173 = 0.1e1 / t175;
t258 = t173 - 0.1e1;
t180 = t187 ^ 2;
t243 = t180 * t182;
t146 = t158 * t243 + 0.1e1;
t234 = qJD(1) * t191;
t210 = t182 * t187 * t234;
t235 = qJD(1) * t190;
t219 = t187 * t235;
t231 = qJD(3) * t191;
t139 = (-(-t186 * t231 - t219) * t177 + t231 * t244) * t173;
t246 = t171 * t190;
t130 = (t139 * t191 - qJD(3)) * t246 + (-t219 + (t139 - t231) * t186) * t170;
t256 = t130 * t157 * t158;
t257 = (-t243 * t256 + (-t180 * t186 * t232 + t210) * t158) / t146 ^ 2;
t150 = t148 * t149;
t255 = t133 * t150;
t132 = qJD(6) * t154 - t142 * t188 + t143 * t184;
t147 = t207 ^ 2;
t138 = t147 * t149 + 0.1e1;
t254 = 0.1e1 / t138 ^ 2 * (-t132 * t251 - t147 * t255);
t253 = t139 * t170;
t252 = t139 * t190;
t239 = t187 * t190;
t164 = t204 * t239;
t250 = t149 * t164;
t249 = t158 * t187;
t248 = t158 * t190;
t162 = t173 * t262;
t247 = t162 * t191;
t245 = t177 * t182;
t242 = t186 * t191;
t236 = qJD(1) * t187;
t233 = qJD(3) * t186;
t226 = 0.2e1 * t256;
t225 = 0.2e1 * t254;
t224 = -0.2e1 * t150 * t207;
t200 = qJD(3) * (-t190 - t263) * t177;
t223 = -0.2e1 * (-t178 * t210 + t183 * t200) / t175 ^ 2;
t222 = t190 * t257;
t221 = t158 * t239;
t220 = t191 * t245;
t218 = t190 * t234;
t216 = t133 * t224;
t215 = t190 * t223;
t212 = t173 * t220;
t211 = t258 * t190 * t170;
t209 = t217 * t187;
t168 = t185 * t242 + t240;
t169 = t186 * t238 - t241;
t206 = t168 * t188 - t169 * t184;
t156 = t168 * t184 + t169 * t188;
t199 = -t187 * t213 + t190 * t231;
t163 = t205 * t239;
t144 = 0.1e1 / t146;
t141 = t189 * t199 - t191 * t203;
t140 = t185 * t199 + t191 * t202;
t136 = 0.1e1 / t138;
t135 = (t171 * t212 + t211) * t187;
t134 = -t170 * t242 - t246 + (t170 * t186 + t171 * t237) * t162;
t131 = t223 * t262 + (-qJD(1) * t209 + 0.2e1 * t191 * t200) * t173;
t1 = [t187 * t177 * t215 + (-qJD(3) * t209 + t177 * t218) * t173, 0, t131, 0, 0, 0; (-0.2e1 * t157 * t222 + (-t157 * t233 + (-qJD(1) * t135 - t130) * t248) * t144) * t191 + (0.2e1 * t158 * t222 * t135 + (-((-t139 * t212 - t233 * t258 + t215) * t170 + (t220 * t223 - t252 + (t252 + (-0.2e1 * t190 - t263) * t231) * t173) * t171) * t221 + (t158 * t233 + t190 * t226) * t135 + (-t157 + ((t180 - t183) * t173 * t171 * t245 - t191 * t211) * t158) * t235) * t144) * t187, 0, 0.2e1 * (t134 * t248 + t157 * t186) * t187 * t257 + ((-t157 * t234 + (qJD(3) * t134 + t130) * t249) * t186 + (-t187 * qJD(3) * t157 - (t131 * t171 * t191 - t170 * t231 - t247 * t253 + t253 + (qJD(3) * t170 - t171 * t236) * t162) * t221 + (-t158 * t234 + t187 * t226) * t134 - ((t131 + t236) * t170 + ((0.1e1 - t247) * qJD(3) + (t162 - t191) * t139) * t171) * t186 * t249) * t190) * t144, 0, 0, 0; (t148 * t206 - t156 * t251) * t225 + ((qJD(6) * t156 - t140 * t188 + t141 * t184) * t148 + t156 * t216 + (t206 * t133 + (qJD(6) * t206 + t140 * t184 + t141 * t188) * t207 - t156 * t132) * t149) * t136, 0 (-t148 * t163 - t207 * t250) * t225 + (-t132 * t250 + (-t149 * t163 + t164 * t224) * t133 + t259 * t218 + (-t259 * t233 + ((t148 * t261 + t251 * t260) * t188 + (t148 * t260 - t251 * t261) * t184) * t190) * t187) * t136 (t148 * t154 + t207 * t251) * t225 + (-t133 * t148 - t207 * t216 + (0.2e1 * t132 * t207 + t133 * t154) * t149) * t136, 0, -0.2e1 * t254 - 0.2e1 * (t132 * t149 * t136 - (-t136 * t255 - t149 * t254) * t207) * t207;];
JaD_rot  = t1;
