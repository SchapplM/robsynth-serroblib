% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR3
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR3_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobiaD_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:24
% EndTime: 2019-02-26 20:05:24
% DurationCPUTime: 0.53s
% Computational Cost: add. (1687->72), mult. (5242->175), div. (282->12), fcn. (6739->15), ass. (0->87)
t200 = sin(pkin(7));
t197 = t200 ^ 2;
t246 = 0.2e1 * t197;
t198 = sin(pkin(13));
t202 = cos(pkin(13));
t206 = sin(qJ(3));
t208 = cos(qJ(3));
t192 = t206 * t198 - t208 * t202;
t204 = cos(pkin(7));
t177 = t192 * t204;
t199 = sin(pkin(12));
t203 = cos(pkin(12));
t209 = cos(qJ(2));
t205 = cos(pkin(6));
t207 = sin(qJ(2));
t234 = t205 * t207;
t218 = t199 * t234 - t203 * t209;
t233 = t205 * t209;
t219 = t199 * t233 + t203 * t207;
t221 = t208 * t198 + t206 * t202;
t201 = sin(pkin(6));
t225 = t199 * t200 * t201;
t158 = t177 * t219 - t192 * t225 + t218 * t221;
t153 = t158 ^ 2;
t178 = t221 * t204;
t216 = -t178 * t219 + t192 * t218 + t221 * t225;
t155 = 0.1e1 / t216 ^ 2;
t245 = t153 * t155;
t220 = -t199 * t207 + t203 * t233;
t236 = t201 * t204;
t171 = t220 * t200 + t203 * t236;
t237 = t200 * t209;
t186 = -t201 * t237 + t205 * t204;
t168 = atan2(t171, t186);
t163 = sin(t168);
t164 = cos(t168);
t151 = t163 * t171 + t164 * t186;
t148 = 0.1e1 / t151;
t154 = 0.1e1 / t216;
t183 = 0.1e1 / t186;
t149 = 0.1e1 / t151 ^ 2;
t184 = 0.1e1 / t186 ^ 2;
t172 = t199 * t236 + t200 * t219;
t170 = t172 ^ 2;
t145 = t170 * t149 + 0.1e1;
t182 = t218 * qJD(2);
t187 = -t199 * t209 - t203 * t234;
t180 = t187 * qJD(2);
t235 = t201 * t207;
t239 = t171 * t184;
t223 = t235 * t239;
t169 = t171 ^ 2;
t167 = t169 * t184 + 0.1e1;
t165 = 0.1e1 / t167;
t240 = t165 * t200;
t138 = (-qJD(2) * t223 + t180 * t183) * t240;
t222 = -t163 * t186 + t164 * t171;
t230 = qJD(2) * t201;
t224 = t207 * t230;
t136 = (t163 * t180 + t164 * t224) * t200 + t222 * t138;
t243 = t136 * t148 * t149;
t244 = (-t172 * t149 * t182 * t200 - t170 * t243) / t145 ^ 2;
t174 = qJD(3) * t177;
t181 = t219 * qJD(2);
t190 = t192 * qJD(3);
t191 = t221 * qJD(3);
t147 = t174 * t219 + t182 * t178 + t181 * t192 - t190 * t225 + t191 * t218;
t156 = t154 * t155;
t242 = t147 * t156;
t162 = t178 * t218 + t192 * t219;
t241 = t158 * t162;
t142 = 0.1e1 + t245;
t175 = t204 * t191;
t146 = t175 * t219 - t182 * t177 + t181 * t221 - t190 * t218 - t191 * t225;
t227 = t158 * t155 * t146;
t228 = 0.2e1 * (-t153 * t242 + t227) / t142 ^ 2;
t226 = t184 * t197 * t207;
t217 = -t183 * t187 + t223;
t185 = t183 * t184;
t179 = t220 * qJD(2);
t161 = t177 * t218 - t219 * t221;
t143 = 0.1e1 / t145;
t140 = 0.1e1 / t142;
t139 = t217 * t240;
t137 = (t163 * t187 + t164 * t235) * t200 - t222 * t139;
t134 = t217 / t167 ^ 2 * (-t169 * t185 * t224 + t180 * t239) * t246 + (-t179 * t183 * t200 + (-t180 * t226 + (-t187 * t226 + (t185 * t201 * t207 ^ 2 * t246 - t184 * t237) * t171) * qJD(2)) * t201) * t165;
t1 = [0, t134, 0, 0, 0, 0; 0 (-(t151 * t139 * t138 + t222 * t134) * t149 * t143 + 0.2e1 * (t143 * t243 + t149 * t244) * t137) * t172 + (0.2e1 * t218 * t148 * t244 + (-t181 * t148 + (t218 * t136 + t137 * t182 + (-(t138 * t187 - t139 * t180 + t209 * t230) * t164 - (-t179 + (qJD(2) * t139 - t138) * t235) * t163) * t172) * t149) * t143) * t200, 0, 0, 0, 0; 0 (-t154 * t161 - t155 * t241) * t228 + ((t175 * t218 + t181 * t177 + t182 * t221 + t190 * t219) * t154 - 0.2e1 * t241 * t242 + (-t161 * t147 + (-t174 * t218 + t181 * t178 - t182 * t192 + t191 * t219) * t158 + t162 * t146) * t155) * t140 (-t154 * t216 - t245) * t228 + (0.2e1 * t227 + (-0.2e1 * t153 * t156 - t155 * t216 + t154) * t147) * t140, 0, 0, 0;];
JaD_rot  = t1;
