% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR8_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:40
% EndTime: 2019-02-26 22:07:40
% DurationCPUTime: 0.64s
% Computational Cost: add. (1349->91), mult. (4347->200), div. (672->14), fcn. (5566->11), ass. (0->94)
t183 = sin(qJ(2));
t184 = sin(qJ(1));
t186 = cos(qJ(2));
t187 = cos(qJ(1));
t236 = cos(pkin(6));
t205 = t187 * t236;
t165 = t183 * t205 + t184 * t186;
t206 = t184 * t236;
t166 = t187 * t183 + t186 * t206;
t145 = t165 * qJD(1) + t166 * qJD(2);
t182 = sin(qJ(3));
t185 = cos(qJ(3));
t181 = sin(pkin(6));
t220 = qJD(1) * t181;
t208 = t187 * t220;
t203 = t183 * t206;
t221 = t187 * t186;
t168 = -t203 + t221;
t225 = t181 * t184;
t155 = t168 * t185 + t182 * t225;
t217 = qJD(3) * t155;
t135 = -t145 * t182 - t185 * t208 + t217;
t195 = -t168 * t182 + t185 * t225;
t148 = 0.1e1 / t195;
t149 = 0.1e1 / t195 ^ 2;
t150 = t148 * t149;
t240 = 0.2e1 * t155 * t150 * t135;
t146 = t166 * qJD(1) + t165 * qJD(2);
t202 = t186 * t205;
t222 = t184 * t183;
t164 = -t202 + t222;
t162 = t164 ^ 2;
t177 = 0.1e1 / t181 ^ 2;
t179 = 0.1e1 / t186 ^ 2;
t161 = t162 * t177 * t179 + 0.1e1;
t178 = 0.1e1 / t186;
t180 = t178 * t179;
t219 = qJD(2) * t183;
t233 = (t146 * t164 * t179 + t162 * t180 * t219) * t177 / t161 ^ 2;
t239 = -0.2e1 * t233;
t176 = 0.1e1 / t181;
t238 = t164 * t176;
t224 = t181 * t186;
t160 = atan2(t164, t224);
t156 = sin(t160);
t157 = cos(t160);
t158 = 0.1e1 / t161;
t210 = t178 * t238;
t237 = (t157 * t210 - t156) * t158 + t156;
t140 = t156 * t164 + t157 * t224;
t137 = 0.1e1 / t140;
t138 = 0.1e1 / t140 ^ 2;
t163 = t166 ^ 2;
t133 = t163 * t138 + 0.1e1;
t201 = qJD(2) * t236 + qJD(1);
t218 = qJD(2) * t186;
t144 = -qJD(1) * t202 - t187 * t218 + t201 * t222;
t232 = t138 * t166;
t207 = t179 * t219;
t194 = (t146 * t178 + t164 * t207) * t176;
t129 = t158 * t194;
t198 = -t156 * t224 + t157 * t164;
t211 = t157 * t181 * t183;
t125 = -qJD(2) * t211 + t198 * t129 + t156 * t146;
t234 = t125 * t137 * t138;
t235 = (-t144 * t232 - t163 * t234) / t133 ^ 2;
t231 = t149 * t155;
t151 = t155 ^ 2;
t230 = t150 * t151;
t229 = t151 * t149;
t228 = t156 * t166;
t227 = t157 * t166;
t226 = t179 * t183;
t223 = t181 * t187;
t216 = -0.2e1 * t235;
t215 = -0.2e1 * t234;
t143 = 0.1e1 + t229;
t136 = t195 * qJD(3) - t145 * t185 + t182 * t208;
t213 = t136 * t231;
t214 = 0.2e1 * (t135 * t230 + t213) / t143 ^ 2;
t209 = t184 * t220;
t204 = t178 * t239;
t199 = -t185 * t148 - t182 * t231;
t197 = t164 * t226 + t165 * t178;
t196 = t165 * t182 + t185 * t223;
t153 = -t165 * t185 + t182 * t223;
t147 = -qJD(1) * t203 - t184 * t219 + t201 * t221;
t141 = 0.1e1 / t143;
t131 = 0.1e1 / t133;
t130 = t197 * t176 * t158;
t128 = t237 * t166;
t126 = t198 * t130 + t156 * t165 - t211;
t124 = (t197 * t239 + (t146 * t226 + t147 * t178 + (t165 * t226 + (0.2e1 * t180 * t183 ^ 2 + t178) * t164) * qJD(2)) * t158) * t176;
t1 = [(t166 * t204 + (-t144 * t178 + t166 * t207) * t158) * t176, t124, 0, 0, 0, 0; t164 * t137 * t216 + (t146 * t137 + (-t125 * t164 - t128 * t144) * t138) * t131 + (t128 * t215 * t131 + (t128 * t216 + ((-t129 * t158 * t210 + 0.2e1 * t233) * t228 + (t204 * t238 + t129 + (-t129 + t194) * t158) * t227 - t237 * t144) * t131) * t138) * t166, 0.2e1 * (-t126 * t232 + t137 * t168) * t235 + (t126 * t166 * t215 + t145 * t137 + (t168 * t125 - t126 * t144 + (-t181 * t218 + t124 * t164 + t130 * t146 + (-t130 * t224 + t165) * t129) * t227 + (-t129 * t130 * t164 + t147 + (-t124 * t186 + (qJD(2) * t130 + t129) * t183) * t181) * t228) * t138) * t131, 0, 0, 0, 0; (t148 * t153 - t196 * t231) * t214 + (-(t196 * qJD(3) - t147 * t185 - t182 * t209) * t148 + t196 * t240 + (-t153 * t135 - (t153 * qJD(3) - t147 * t182 + t185 * t209) * t155 + t196 * t136) * t149) * t141, t199 * t166 * t214 + (t199 * t144 + ((-qJD(3) * t148 + t240) * t182 + (t136 * t182 + (t135 + t217) * t185) * t149) * t166) * t141 (t148 * t195 + t229) * t214 + (-0.2e1 * t213 + (-t149 * t195 + t148 - 0.2e1 * t230) * t135) * t141, 0, 0, 0;];
JaD_rot  = t1;
