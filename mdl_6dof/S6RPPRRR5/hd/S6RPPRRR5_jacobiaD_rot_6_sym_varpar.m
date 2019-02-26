% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:37:04
% EndTime: 2019-02-26 20:37:05
% DurationCPUTime: 0.71s
% Computational Cost: add. (3075->96), mult. (3810->207), div. (753->12), fcn. (4455->9), ass. (0->97)
t178 = sin(qJ(1));
t176 = qJ(4) + qJ(5);
t171 = sin(t176);
t167 = 0.1e1 / t171 ^ 2;
t172 = cos(t176);
t170 = t172 ^ 2;
t224 = t167 * t170;
t200 = 0.1e1 + t224;
t240 = t178 * t200;
t174 = t178 ^ 2;
t163 = t174 * t224 + 0.1e1;
t159 = 0.1e1 / t163;
t166 = 0.1e1 / t171;
t180 = cos(qJ(1));
t210 = qJD(1) * t180;
t201 = t172 * t210;
t173 = qJD(4) + qJD(5);
t218 = t173 * t178;
t204 = t167 * t218;
t134 = ((-t171 * t218 + t201) * t166 - t170 * t204) * t159;
t239 = -t134 - t218;
t215 = t178 * t172;
t161 = atan2(t215, t171);
t158 = cos(t161);
t157 = sin(t161);
t205 = t157 * t215;
t144 = t158 * t171 + t205;
t141 = 0.1e1 / t144;
t179 = cos(qJ(6));
t212 = t179 * t180;
t203 = t171 * t212;
t177 = sin(qJ(6));
t214 = t178 * t177;
t156 = t203 - t214;
t150 = 0.1e1 / t156;
t142 = 0.1e1 / t144 ^ 2;
t151 = 0.1e1 / t156 ^ 2;
t238 = t159 - 0.1e1;
t175 = t180 ^ 2;
t223 = t170 * t175;
t137 = t142 * t223 + 0.1e1;
t222 = t170 * t178;
t193 = t210 * t222;
t220 = t172 * t173;
t226 = t158 * t172;
t129 = (t134 * t178 + t173) * t226 + (t239 * t171 + t201) * t157;
t236 = t129 * t141 * t142;
t237 = (-t223 * t236 + (-t171 * t175 * t220 - t193) * t142) / t137 ^ 2;
t195 = qJD(1) * t171 + qJD(6);
t217 = t173 * t180;
t202 = t172 * t217;
t138 = -qJD(6) * t203 - t177 * t202 - t179 * t210 + t195 * t214;
t213 = t178 * t179;
t216 = t177 * t180;
t155 = t171 * t216 + t213;
t149 = t155 ^ 2;
t148 = t149 * t151 + 0.1e1;
t229 = t151 * t155;
t196 = qJD(6) * t171 + qJD(1);
t139 = -t196 * t216 + (-t195 * t178 + t202) * t179;
t234 = t139 * t150 * t151;
t235 = (-t138 * t229 - t149 * t234) / t148 ^ 2;
t169 = t172 * t170;
t225 = t166 * t172;
t191 = t166 * t167 * t169 + t225;
t233 = (-t191 * t174 * t173 + t167 * t193) / t163 ^ 2;
t232 = t142 * t172;
t231 = t142 * t180;
t230 = t150 * t177;
t228 = t155 * t179;
t227 = t157 * t178;
t221 = t171 * t173;
t219 = t172 * t180;
t211 = qJD(1) * t178;
t209 = -0.2e1 * t236;
t208 = 0.2e1 * t235;
t207 = 0.2e1 * t233;
t206 = t142 * t219;
t199 = -0.2e1 * t172 * t237;
t198 = 0.2e1 * t155 * t234;
t197 = -0.2e1 * t166 * t233;
t194 = t158 * t159 * t166 * t170;
t192 = t200 * t180;
t190 = t195 * t180;
t189 = t151 * t228 - t230;
t188 = t189 * t180;
t154 = -t171 * t213 - t216;
t153 = -t171 * t214 + t212;
t146 = 0.1e1 / t148;
t145 = t159 * t240;
t135 = 0.1e1 / t137;
t133 = (-t238 * t172 * t157 + t178 * t194) * t180;
t131 = -t171 * t227 + t226 - (-t157 * t171 + t158 * t215) * t145;
t130 = t207 * t240 + (-qJD(1) * t192 + 0.2e1 * t191 * t218) * t159;
t127 = t172 * t188 * t208 + (t188 * t221 + (t189 * t211 + ((qJD(6) * t150 + t198) * t179 + (t138 * t179 + (qJD(6) * t155 - t139) * t177) * t151) * t180) * t172) * t146;
t126 = 0.2e1 * (-t131 * t232 - t141 * t171) * t180 * t237 + ((-t141 * t211 + (-t131 * t173 - t129) * t231) * t171 + (t141 * t217 + (t130 * t158 * t178 + t239 * t157 - (-t134 * t227 - t157 * t173 + t158 * t210) * t145) * t206 + (-t142 * t211 + t180 * t209) * t131 + ((-t130 - t210) * t157 + ((t145 * t178 - 0.1e1) * t173 + (t145 - t178) * t134) * t158) * t171 * t231) * t172) * t135;
t1 = [t197 * t219 + (-t173 * t192 - t211 * t225) * t159, 0, 0, t130, t130, 0; (t141 * t199 + (-t141 * t221 + (-qJD(1) * t133 - t129) * t232) * t135) * t178 + (t142 * t199 * t133 + (((t172 * t207 - t221 + (-t134 * t166 * t222 + t221) * t159) * t157 + (t197 * t222 + t134 * t172 + (-t169 * t204 + (-t134 - 0.2e1 * t218) * t172) * t159) * t158) * t206 + (-t142 * t221 + t172 * t209) * t133 + (t141 + ((-t174 + t175) * t194 + t238 * t205) * t142) * t172 * qJD(1)) * t135) * t180, 0, 0, t126, t126, 0; (-t150 * t153 + t154 * t229) * t208 + (t154 * t198 - t196 * t150 * t213 + (-t173 * t215 - t190) * t230 + (t154 * t138 - t153 * t139 + t190 * t228 - (t196 * t177 - t179 * t220) * t155 * t178) * t151) * t146, 0, 0, t127, t127, -0.2e1 * t235 + 0.2e1 * (-t138 * t146 * t151 + (-t146 * t234 - t151 * t235) * t155) * t155;];
JaD_rot  = t1;
