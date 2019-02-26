% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:31:16
% EndTime: 2019-02-26 22:31:17
% DurationCPUTime: 0.83s
% Computational Cost: add. (8394->97), mult. (5101->205), div. (1026->12), fcn. (5942->9), ass. (0->96)
t182 = sin(qJ(1));
t243 = 0.2e1 * t182;
t180 = t182 ^ 2;
t178 = qJ(2) + qJ(3) + qJ(4);
t173 = sin(t178);
t169 = t173 ^ 2;
t174 = cos(t178);
t171 = 0.1e1 / t174 ^ 2;
t228 = t169 * t171;
t164 = t180 * t228 + 0.1e1;
t162 = 0.1e1 / t164;
t170 = 0.1e1 / t174;
t183 = cos(qJ(1));
t214 = qJD(1) * t183;
t204 = t173 * t214;
t177 = qJD(2) + qJD(3) + qJD(4);
t220 = t177 * t182;
t207 = t171 * t220;
t136 = (-(-t174 * t220 - t204) * t170 + t169 * t207) * t162;
t242 = t136 - t220;
t179 = pkin(11) + qJ(6);
t175 = sin(t179);
t217 = t182 * t175;
t176 = cos(t179);
t221 = t176 * t183;
t158 = t174 * t221 + t217;
t218 = t182 * t173;
t161 = atan2(-t218, -t174);
t160 = cos(t161);
t159 = sin(t161);
t208 = t159 * t218;
t146 = -t160 * t174 - t208;
t143 = 0.1e1 / t146;
t152 = 0.1e1 / t158;
t144 = 0.1e1 / t146 ^ 2;
t153 = 0.1e1 / t158 ^ 2;
t241 = t162 - 0.1e1;
t230 = t160 * t173;
t131 = (-t136 * t182 + t177) * t230 + (t242 * t174 - t204) * t159;
t240 = t131 * t143 * t144;
t193 = t174 * t217 + t221;
t219 = t177 * t183;
t206 = t173 * t219;
t137 = t193 * qJD(1) - t158 * qJD(6) + t175 * t206;
t216 = t182 * t176;
t222 = t175 * t183;
t157 = t174 * t222 - t216;
t151 = t157 ^ 2;
t150 = t151 * t153 + 0.1e1;
t233 = t153 * t157;
t198 = -qJD(1) * t174 + qJD(6);
t199 = qJD(6) * t174 - qJD(1);
t138 = -t199 * t222 + (t198 * t182 - t206) * t176;
t238 = t138 * t152 * t153;
t239 = (-t137 * t233 - t151 * t238) / t150 ^ 2;
t168 = t173 * t169;
t225 = t170 * t173;
t192 = t177 * (t168 * t170 * t171 + t225);
t226 = t169 * t182;
t196 = t214 * t226;
t237 = (t171 * t196 + t180 * t192) / t164 ^ 2;
t236 = t144 * t173;
t235 = t144 * t183;
t234 = t152 * t175;
t232 = t157 * t176;
t231 = t159 * t182;
t229 = t169 * t170;
t181 = t183 ^ 2;
t227 = t169 * t181;
t224 = t173 * t183;
t223 = t174 * t177;
t215 = qJD(1) * t182;
t141 = t144 * t227 + 0.1e1;
t213 = 0.2e1 * (-t227 * t240 + (t173 * t181 * t223 - t196) * t144) / t141 ^ 2;
t212 = 0.2e1 * t240;
t211 = -0.2e1 * t239;
t210 = t157 * t238;
t209 = t144 * t224;
t203 = 0.1e1 + t228;
t202 = t173 * t213;
t201 = -0.2e1 * t173 * t237;
t200 = t237 * t243;
t197 = t160 * t162 * t229;
t195 = t203 * t183;
t194 = t153 * t232 - t234;
t191 = t177 * t218 + t198 * t183;
t156 = -t174 * t216 + t222;
t148 = 0.1e1 / t150;
t147 = t203 * t182 * t162;
t139 = 0.1e1 / t141;
t135 = (t241 * t173 * t159 - t182 * t197) * t183;
t134 = -t174 * t231 + t230 + (t159 * t174 - t160 * t218) * t147;
t132 = -t203 * t200 + (qJD(1) * t195 + t192 * t243) * t162;
t129 = t194 * t211 * t224 + (t194 * t174 * t219 + (-t194 * t215 + ((-qJD(6) * t152 - 0.2e1 * t210) * t176 + (-t137 * t176 + (-qJD(6) * t157 + t138) * t175) * t153) * t183) * t173) * t148;
t128 = (t134 * t236 - t143 * t174) * t183 * t213 + ((-t143 * t215 + (-t134 * t177 - t131) * t235) * t174 + (-t143 * t219 - (-t132 * t160 * t182 - t242 * t159 + (t136 * t231 - t159 * t177 - t160 * t214) * t147) * t209 + (t144 * t215 + t183 * t212) * t134 - ((t132 - t214) * t159 + ((-t147 * t182 + 0.1e1) * t177 + (t147 - t182) * t136) * t160) * t174 * t235) * t173) * t139;
t1 = [t170 * t183 * t201 + (t177 * t195 - t215 * t225) * t162, t132, t132, t132, 0, 0; (t143 * t202 + (-t143 * t223 + (qJD(1) * t135 + t131) * t236) * t139) * t182 + (t144 * t202 * t135 + (-((t201 - t223 + (t136 * t170 * t226 + t223) * t162) * t159 + (t200 * t229 - t136 * t173 + (-t168 * t207 + (t136 - 0.2e1 * t220) * t173) * t162) * t160) * t209 + (-t144 * t223 + t173 * t212) * t135 + (-t143 + ((-t180 + t181) * t197 + t241 * t208) * t144) * t173 * qJD(1)) * t139) * t183, t128, t128, t128, 0, 0; 0.2e1 * (t152 * t193 + t156 * t233) * t239 + (0.2e1 * t156 * t210 - t199 * t152 * t216 + t191 * t234 + (-t199 * t157 * t217 + t156 * t137 + t138 * t193 - t191 * t232) * t153) * t148, t129, t129, t129, 0, t211 + 0.2e1 * (-t137 * t148 * t153 + (-t148 * t238 - t153 * t239) * t157) * t157;];
JaD_rot  = t1;
