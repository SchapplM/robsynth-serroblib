% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:31:42
% EndTime: 2019-02-26 22:31:43
% DurationCPUTime: 0.74s
% Computational Cost: add. (7389->95), mult. (5101->207), div. (1026->12), fcn. (5942->9), ass. (0->95)
t180 = sin(qJ(1));
t176 = qJ(2) + qJ(3) + qJ(4);
t173 = sin(t176);
t169 = 0.1e1 / t173 ^ 2;
t174 = cos(t176);
t172 = t174 ^ 2;
t226 = t169 * t172;
t201 = 0.1e1 + t226;
t241 = t180 * t201;
t177 = t180 ^ 2;
t166 = t177 * t226 + 0.1e1;
t164 = 0.1e1 / t166;
t168 = 0.1e1 / t173;
t182 = cos(qJ(1));
t213 = qJD(1) * t182;
t202 = t174 * t213;
t175 = qJD(2) + qJD(3) + qJD(4);
t221 = t175 * t180;
t204 = t169 * t221;
t138 = ((t173 * t221 - t202) * t168 + t172 * t204) * t164;
t240 = -t138 + t221;
t197 = qJD(1) * t173 + qJD(6);
t220 = t175 * t182;
t239 = -t174 * t220 + t197 * t180;
t217 = t180 * t174;
t163 = atan2(-t217, t173);
t158 = cos(t163);
t157 = sin(t163);
t206 = t157 * t217;
t148 = t158 * t173 - t206;
t145 = 0.1e1 / t148;
t181 = cos(qJ(6));
t216 = t180 * t181;
t179 = sin(qJ(6));
t218 = t179 * t182;
t160 = t173 * t218 + t216;
t154 = 0.1e1 / t160;
t146 = 0.1e1 / t148 ^ 2;
t155 = 0.1e1 / t160 ^ 2;
t238 = t164 - 0.1e1;
t229 = t158 * t174;
t133 = (-t138 * t180 + t175) * t229 + (t240 * t173 - t202) * t157;
t237 = t133 * t145 * t146;
t198 = qJD(6) * t173 + qJD(1);
t193 = t198 * t182;
t143 = t179 * t193 + t239 * t181;
t215 = t181 * t182;
t219 = t179 * t180;
t159 = -t173 * t215 + t219;
t153 = t159 ^ 2;
t152 = t153 * t155 + 0.1e1;
t231 = t155 * t159;
t144 = -t239 * t179 + t181 * t193;
t234 = t144 * t154 * t155;
t236 = (t143 * t231 - t153 * t234) / t152 ^ 2;
t171 = t174 * t172;
t227 = t168 * t174;
t191 = t175 * (-t168 * t169 * t171 - t227);
t224 = t172 * t180;
t195 = t213 * t224;
t235 = (t169 * t195 + t177 * t191) / t166 ^ 2;
t233 = t146 * t174;
t232 = t146 * t182;
t230 = t157 * t180;
t228 = t159 * t179;
t178 = t182 ^ 2;
t225 = t172 * t178;
t223 = t173 * t175;
t222 = t174 * t175;
t214 = qJD(1) * t180;
t141 = t146 * t225 + 0.1e1;
t212 = 0.2e1 * (-t225 * t237 + (-t173 * t178 * t222 - t195) * t146) / t141 ^ 2;
t211 = 0.2e1 * t237;
t210 = 0.2e1 * t236;
t209 = -0.2e1 * t235;
t208 = t174 * t235;
t207 = t174 * t232;
t205 = t168 * t224;
t200 = t174 * t212;
t199 = 0.2e1 * t159 * t234;
t196 = t158 * t164 * t168 * t172;
t194 = t201 * t182;
t192 = t154 * t181 + t155 * t228;
t190 = t192 * t182;
t162 = -t173 * t219 + t215;
t161 = t173 * t216 + t218;
t150 = 0.1e1 / t152;
t149 = t164 * t241;
t139 = 0.1e1 / t141;
t137 = (t238 * t174 * t157 + t180 * t196) * t182;
t135 = t173 * t230 + t229 + (-t157 * t173 - t158 * t217) * t149;
t134 = t209 * t241 + (qJD(1) * t194 + 0.2e1 * t180 * t191) * t164;
t131 = t174 * t190 * t210 + (t190 * t223 + (t192 * t214 + ((qJD(6) * t154 + t199) * t179 + (-t143 * t179 + (-qJD(6) * t159 + t144) * t181) * t155) * t182) * t174) * t150;
t130 = (t135 * t233 + t145 * t173) * t182 * t212 + ((t145 * t214 + (t135 * t175 + t133) * t232) * t173 + (-t145 * t220 - (-t134 * t158 * t180 + t240 * t157 + (t138 * t230 - t157 * t175 - t158 * t213) * t149) * t207 + (t146 * t214 + t182 * t211) * t135 - ((-t134 + t213) * t157 + ((t149 * t180 - 0.1e1) * t175 + (-t149 + t180) * t138) * t158) * t173 * t232) * t174) * t139;
t1 = [0.2e1 * t168 * t182 * t208 + (t175 * t194 + t214 * t227) * t164, t134, t134, t134, 0, 0; (t145 * t200 + (t145 * t223 + (qJD(1) * t137 + t133) * t233) * t139) * t180 + (t146 * t200 * t137 + (-((-0.2e1 * t208 + t223 + (-t138 * t205 - t223) * t164) * t157 + (t205 * t209 - t138 * t174 + (-t171 * t204 + (t138 - 0.2e1 * t221) * t174) * t164) * t158) * t207 + (t146 * t223 + t174 * t211) * t137 + (-t145 + ((t177 - t178) * t196 + t238 * t206) * t146) * t174 * qJD(1)) * t139) * t182, t130, t130, t130, 0, 0; (-t154 * t161 + t162 * t231) * t210 + (t162 * t199 + (-t162 * t143 - t161 * t144 + t198 * t159 * t216 - (-t175 * t217 - t197 * t182) * t228) * t155 + (t197 * t215 + (-t198 * t179 + t181 * t222) * t180) * t154) * t150, t131, t131, t131, 0, -0.2e1 * t236 + 0.2e1 * (t143 * t150 * t155 + (-t150 * t234 - t155 * t236) * t159) * t159;];
JaD_rot  = t1;
