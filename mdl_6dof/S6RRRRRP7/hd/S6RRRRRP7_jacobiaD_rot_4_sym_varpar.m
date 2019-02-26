% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRP7_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:43:03
% EndTime: 2019-02-26 22:43:04
% DurationCPUTime: 0.66s
% Computational Cost: add. (2050->93), mult. (4665->202), div. (686->14), fcn. (5929->11), ass. (0->94)
t202 = cos(qJ(2));
t203 = cos(qJ(1));
t252 = cos(pkin(6));
t224 = t203 * t252;
t222 = t202 * t224;
t200 = sin(qJ(2));
t201 = sin(qJ(1));
t238 = t200 * t201;
t179 = -t222 + t238;
t199 = sin(pkin(6));
t240 = t199 * t202;
t173 = atan2(-t179, -t240);
t171 = sin(t173);
t172 = cos(t173);
t177 = t179 ^ 2;
t193 = 0.1e1 / t199 ^ 2;
t196 = 0.1e1 / t202 ^ 2;
t176 = t177 * t193 * t196 + 0.1e1;
t174 = 0.1e1 / t176;
t192 = 0.1e1 / t199;
t195 = 0.1e1 / t202;
t227 = t179 * t192 * t195;
t253 = (t172 * t227 - t171) * t174 + t171;
t155 = -t171 * t179 - t172 * t240;
t152 = 0.1e1 / t155;
t225 = t201 * t252;
t223 = t200 * t225;
t237 = t203 * t202;
t183 = -t223 + t237;
t198 = qJ(3) + qJ(4);
t190 = sin(t198);
t191 = cos(t198);
t241 = t199 * t201;
t170 = t183 * t191 + t190 * t241;
t164 = 0.1e1 / t170;
t153 = 0.1e1 / t155 ^ 2;
t165 = 0.1e1 / t170 ^ 2;
t210 = -t200 * t224 - t201 * t202;
t211 = -t203 * t200 - t202 * t225;
t161 = -t211 * qJD(1) - t210 * qJD(2);
t235 = qJD(2) * t200;
t226 = t196 * t235;
t212 = t161 * t195 + t179 * t226;
t243 = t174 * t192;
t144 = t212 * t243;
t216 = t171 * t240 - t172 * t179;
t228 = t172 * t199 * t200;
t140 = qJD(2) * t228 + t216 * t144 - t171 * t161;
t251 = t140 * t152 * t153;
t194 = qJD(3) + qJD(4);
t236 = qJD(1) * t199;
t213 = -t183 * t194 + t203 * t236;
t160 = t210 * qJD(1) + t211 * qJD(2);
t221 = t194 * t241 + t160;
t149 = t213 * t190 + t221 * t191;
t250 = t149 * t164 * t165;
t197 = t195 * t196;
t249 = (t161 * t179 * t196 + t177 * t197 * t235) * t193 / t176 ^ 2;
t248 = t153 * t211;
t148 = t221 * t190 - t213 * t191;
t169 = t183 * t190 - t191 * t241;
t163 = t169 ^ 2;
t158 = t163 * t165 + 0.1e1;
t246 = t165 * t169;
t247 = 0.1e1 / t158 ^ 2 * (t148 * t246 - t163 * t250);
t245 = t171 * t211;
t244 = t172 * t211;
t242 = t196 * t200;
t239 = t199 * t203;
t234 = qJD(2) * t202;
t178 = t211 ^ 2;
t150 = t153 * t178 + 0.1e1;
t219 = qJD(2) * t252 + qJD(1);
t159 = -qJD(1) * t222 - t203 * t234 + t219 * t238;
t233 = 0.2e1 * (t159 * t248 - t178 * t251) / t150 ^ 2;
t232 = 0.2e1 * t251;
t231 = -0.2e1 * t249;
t230 = 0.2e1 * t247;
t229 = t169 * t250;
t162 = -qJD(1) * t223 - t201 * t235 + t219 * t237;
t220 = t194 * t239 - t162;
t217 = t164 * t190 - t191 * t246;
t215 = t179 * t242 - t195 * t210;
t214 = t194 * t210 + t201 * t236;
t168 = t190 * t239 + t191 * t210;
t167 = t190 * t210 - t191 * t239;
t156 = 0.1e1 / t158;
t146 = 0.1e1 / t150;
t145 = t215 * t243;
t143 = t253 * t211;
t141 = t216 * t145 + t171 * t210 + t228;
t139 = (t215 * t231 + (t161 * t242 + t162 * t195 + (-t210 * t242 + (0.2e1 * t197 * t200 ^ 2 + t195) * t179) * qJD(2)) * t174) * t192;
t137 = -0.2e1 * t247 + 0.2e1 * (t148 * t165 * t156 + (-t156 * t250 - t165 * t247) * t169) * t169;
t1 = [(-t211 * t195 * t231 + (-t159 * t195 - t211 * t226) * t174) * t192, t139, 0, 0, 0, 0; t179 * t152 * t233 + (-t161 * t152 + (t140 * t179 + t143 * t159) * t153) * t146 - (t143 * t232 * t146 + (t143 * t233 + ((t144 * t174 * t227 + t231) * t245 + (0.2e1 * t227 * t249 - t144 + (-t212 * t192 + t144) * t174) * t244 - t253 * t159) * t146) * t153) * t211 (-t141 * t248 - t152 * t183) * t233 + (-t141 * t211 * t232 + t160 * t152 + (-t183 * t140 + t141 * t159 + (t199 * t234 - t139 * t179 - t145 * t161 + (t145 * t240 + t210) * t144) * t244 + (t144 * t145 * t179 - t162 + (t139 * t202 + (-qJD(2) * t145 - t144) * t200) * t199) * t245) * t153) * t146, 0, 0, 0, 0; (-t164 * t167 + t168 * t246) * t230 + ((t220 * t190 + t214 * t191) * t164 + 0.2e1 * t168 * t229 + (-t167 * t149 - (-t214 * t190 + t220 * t191) * t169 - t168 * t148) * t165) * t156, -t217 * t211 * t230 + (t217 * t159 - ((-t164 * t194 - 0.2e1 * t229) * t191 + (t148 * t191 + (-t169 * t194 + t149) * t190) * t165) * t211) * t156, t137, t137, 0, 0;];
JaD_rot  = t1;
