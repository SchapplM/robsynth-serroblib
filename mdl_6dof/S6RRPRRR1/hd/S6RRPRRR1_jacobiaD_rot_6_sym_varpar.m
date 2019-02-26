% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:54:01
% EndTime: 2019-02-26 21:54:02
% DurationCPUTime: 0.82s
% Computational Cost: add. (10386->97), mult. (5101->208), div. (1026->12), fcn. (5942->9), ass. (0->95)
t182 = sin(qJ(1));
t244 = 0.2e1 * t182;
t179 = t182 ^ 2;
t176 = qJ(2) + pkin(11) + qJ(4) + qJ(5);
t173 = sin(t176);
t169 = t173 ^ 2;
t174 = cos(t176);
t171 = 0.1e1 / t174 ^ 2;
t229 = t169 * t171;
t166 = t179 * t229 + 0.1e1;
t160 = 0.1e1 / t166;
t170 = 0.1e1 / t174;
t184 = cos(qJ(1));
t215 = qJD(1) * t184;
t205 = t173 * t215;
t178 = qJD(2) + qJD(4) + qJD(5);
t223 = t178 * t182;
t208 = t171 * t223;
t138 = (-(-t174 * t223 - t205) * t170 + t169 * t208) * t160;
t243 = t138 - t223;
t183 = cos(qJ(6));
t217 = t183 * t184;
t181 = sin(qJ(6));
t219 = t182 * t181;
t165 = t174 * t217 + t219;
t220 = t182 * t173;
t155 = atan2(-t220, -t174);
t154 = cos(t155);
t153 = sin(t155);
t209 = t153 * t220;
t146 = -t154 * t174 - t209;
t143 = 0.1e1 / t146;
t157 = 0.1e1 / t165;
t144 = 0.1e1 / t146 ^ 2;
t158 = 0.1e1 / t165 ^ 2;
t242 = t160 - 0.1e1;
t234 = t154 * t173;
t133 = (-t138 * t182 + t178) * t234 + (t243 * t174 - t205) * t153;
t241 = t133 * t143 * t144;
t193 = t174 * t219 + t217;
t222 = t178 * t184;
t207 = t173 * t222;
t147 = t193 * qJD(1) - t165 * qJD(6) + t181 * t207;
t218 = t182 * t183;
t221 = t181 * t184;
t164 = t174 * t221 - t218;
t156 = t164 ^ 2;
t152 = t156 * t158 + 0.1e1;
t232 = t158 * t164;
t199 = -qJD(1) * t174 + qJD(6);
t200 = qJD(6) * t174 - qJD(1);
t148 = -t200 * t221 + (t199 * t182 - t207) * t183;
t236 = t148 * t157 * t158;
t240 = (-t147 * t232 - t156 * t236) / t152 ^ 2;
t168 = t173 * t169;
t226 = t170 * t173;
t192 = t178 * (t168 * t170 * t171 + t226);
t227 = t169 * t182;
t197 = t215 * t227;
t239 = (t171 * t197 + t179 * t192) / t166 ^ 2;
t238 = t144 * t173;
t237 = t144 * t184;
t235 = t153 * t182;
t233 = t157 * t181;
t231 = t164 * t183;
t230 = t169 * t170;
t180 = t184 ^ 2;
t228 = t169 * t180;
t225 = t173 * t184;
t224 = t174 * t178;
t216 = qJD(1) * t182;
t141 = t144 * t228 + 0.1e1;
t214 = 0.2e1 * (-t228 * t241 + (t173 * t180 * t224 - t197) * t144) / t141 ^ 2;
t213 = 0.2e1 * t241;
t212 = -0.2e1 * t240;
t211 = t144 * t225;
t210 = t164 * t236;
t204 = 0.1e1 + t229;
t203 = t173 * t214;
t202 = -0.2e1 * t173 * t239;
t201 = t239 * t244;
t198 = t154 * t160 * t230;
t196 = t204 * t184;
t195 = t199 * t184;
t194 = t158 * t231 - t233;
t163 = -t174 * t218 + t221;
t150 = 0.1e1 / t152;
t149 = t204 * t182 * t160;
t139 = 0.1e1 / t141;
t137 = (t242 * t173 * t153 - t182 * t198) * t184;
t135 = -t174 * t235 + t234 + (t153 * t174 - t154 * t220) * t149;
t134 = -t204 * t201 + (qJD(1) * t196 + t192 * t244) * t160;
t131 = t194 * t212 * t225 + (t194 * t174 * t222 + (-t194 * t216 + ((-qJD(6) * t157 - 0.2e1 * t210) * t183 + (-t147 * t183 + (-qJD(6) * t164 + t148) * t181) * t158) * t184) * t173) * t150;
t130 = (t135 * t238 - t143 * t174) * t184 * t214 + ((-t143 * t216 + (-t135 * t178 - t133) * t237) * t174 + (-t143 * t222 - (-t134 * t154 * t182 - t243 * t153 + (t138 * t235 - t153 * t178 - t154 * t215) * t149) * t211 + (t144 * t216 + t184 * t213) * t135 - ((t134 - t215) * t153 + ((-t149 * t182 + 0.1e1) * t178 + (t149 - t182) * t138) * t154) * t174 * t237) * t173) * t139;
t1 = [t170 * t184 * t202 + (t178 * t196 - t216 * t226) * t160, t134, 0, t134, t134, 0; (t143 * t203 + (-t143 * t224 + (qJD(1) * t137 + t133) * t238) * t139) * t182 + (t144 * t203 * t137 + (-((t202 - t224 + (t138 * t170 * t227 + t224) * t160) * t153 + (t201 * t230 - t138 * t173 + (-t168 * t208 + (t138 - 0.2e1 * t223) * t173) * t160) * t154) * t211 + (-t144 * t224 + t173 * t213) * t137 + (-t143 + ((-t179 + t180) * t198 + t242 * t209) * t144) * t173 * qJD(1)) * t139) * t184, t130, 0, t130, t130, 0; 0.2e1 * (t157 * t193 + t163 * t232) * t240 + (0.2e1 * t163 * t210 - t200 * t157 * t218 + (t178 * t220 + t195) * t233 + (t163 * t147 + t193 * t148 - t195 * t231 - (t173 * t178 * t183 + t200 * t181) * t164 * t182) * t158) * t150, t131, 0, t131, t131, t212 + 0.2e1 * (-t147 * t150 * t158 + (-t150 * t236 - t158 * t240) * t164) * t164;];
JaD_rot  = t1;
