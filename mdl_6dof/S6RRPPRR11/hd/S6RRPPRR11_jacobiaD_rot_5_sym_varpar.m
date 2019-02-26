% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR11_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:19
% EndTime: 2019-02-26 21:34:20
% DurationCPUTime: 0.66s
% Computational Cost: add. (1643->91), mult. (4303->199), div. (668->14), fcn. (5516->11), ass. (0->91)
t181 = sin(qJ(2));
t182 = sin(qJ(1));
t183 = cos(qJ(2));
t184 = cos(qJ(1));
t230 = cos(pkin(6));
t200 = t184 * t230;
t162 = t181 * t200 + t182 * t183;
t180 = sin(pkin(6));
t219 = t180 * t181;
t156 = atan2(-t162, t219);
t152 = sin(t156);
t153 = cos(t156);
t159 = t162 ^ 2;
t175 = 0.1e1 / t180 ^ 2;
t178 = 0.1e1 / t181 ^ 2;
t157 = t159 * t175 * t178 + 0.1e1;
t154 = 0.1e1 / t157;
t174 = 0.1e1 / t180;
t177 = 0.1e1 / t181;
t205 = t162 * t174 * t177;
t231 = (t153 * t205 + t152) * t154 - t152;
t136 = -t152 * t162 + t153 * t219;
t133 = 0.1e1 / t136;
t201 = t182 * t230;
t164 = t184 * t181 + t183 * t201;
t176 = pkin(11) + qJ(5);
t172 = sin(t176);
t173 = cos(t176);
t218 = t180 * t182;
t149 = t164 * t172 + t173 * t218;
t145 = 0.1e1 / t149;
t134 = 0.1e1 / t136 ^ 2;
t146 = 0.1e1 / t149 ^ 2;
t196 = qJD(2) * t230 + qJD(1);
t198 = t181 * t201;
t213 = qJD(2) * t181;
t215 = t184 * t183;
t143 = -qJD(1) * t198 - t182 * t213 + t196 * t215;
t212 = qJD(2) * t183;
t202 = t178 * t212;
t191 = -t143 * t177 + t162 * t202;
t221 = t154 * t174;
t125 = t191 * t221;
t193 = -t152 * t219 - t153 * t162;
t206 = t153 * t180 * t183;
t121 = qJD(2) * t206 + t193 * t125 - t152 * t143;
t229 = t121 * t133 * t134;
t197 = t183 * t200;
t216 = t182 * t181;
t140 = -qJD(1) * t197 - t184 * t212 + t196 * t216;
t214 = qJD(1) * t180;
t203 = t184 * t214;
t130 = t149 * qJD(5) + t140 * t173 + t172 * t203;
t148 = -t164 * t173 + t172 * t218;
t144 = t148 ^ 2;
t139 = t144 * t146 + 0.1e1;
t224 = t146 * t148;
t211 = qJD(5) * t148;
t131 = -t140 * t172 + t173 * t203 - t211;
t227 = t131 * t145 * t146;
t228 = (t130 * t224 - t144 * t227) / t139 ^ 2;
t179 = t177 * t178;
t226 = (t143 * t162 * t178 - t159 * t179 * t212) * t175 / t157 ^ 2;
t165 = -t198 + t215;
t225 = t134 * t165;
t223 = t152 * t165;
t222 = t153 * t165;
t220 = t178 * t183;
t217 = t180 * t184;
t160 = t165 ^ 2;
t129 = t160 * t134 + 0.1e1;
t141 = t162 * qJD(1) + t164 * qJD(2);
t210 = 0.2e1 * (-t141 * t225 - t160 * t229) / t129 ^ 2;
t209 = 0.2e1 * t229;
t208 = 0.2e1 * t228;
t207 = -0.2e1 * t226;
t204 = t182 * t214;
t199 = 0.2e1 * t148 * t227;
t194 = t173 * t145 + t172 * t224;
t161 = -t197 + t216;
t192 = t161 * t177 + t162 * t220;
t151 = -t161 * t172 + t173 * t217;
t150 = t161 * t173 + t172 * t217;
t142 = t164 * qJD(1) + t162 * qJD(2);
t137 = 0.1e1 / t139;
t127 = 0.1e1 / t129;
t126 = t192 * t221;
t124 = t231 * t165;
t122 = t193 * t126 + t152 * t161 + t206;
t120 = (t192 * t207 + (t143 * t220 + t142 * t177 + (-t161 * t220 + (-0.2e1 * t179 * t183 ^ 2 - t177) * t162) * qJD(2)) * t154) * t174;
t1 = [(0.2e1 * t165 * t177 * t226 + (t141 * t177 + t165 * t202) * t154) * t174, t120, 0, 0, 0, 0; t162 * t133 * t210 + (-t143 * t133 + (t121 * t162 + t124 * t141) * t134) * t127 + (t124 * t209 * t127 + (t124 * t210 + (-(-t125 * t154 * t205 + t207) * t223 - (t205 * t207 - t125 + (-t191 * t174 + t125) * t154) * t222 + t231 * t141) * t127) * t134) * t165 (t122 * t225 + t133 * t164) * t210 + (t122 * t165 * t209 + t140 * t133 + (t164 * t121 + t122 * t141 - (-t180 * t213 - t120 * t162 - t126 * t143 + (-t126 * t219 + t161) * t125) * t222 - (t125 * t126 * t162 + t142 + (-t120 * t181 + (-qJD(2) * t126 - t125) * t183) * t180) * t223) * t134) * t127, 0, 0, 0, 0; (-t145 * t150 + t151 * t224) * t208 + ((t151 * qJD(5) + t142 * t173 - t172 * t204) * t145 + t151 * t199 + (-t150 * t131 - (-t150 * qJD(5) - t142 * t172 - t173 * t204) * t148 - t151 * t130) * t146) * t137, t194 * t165 * t208 + (t194 * t141 + ((qJD(5) * t145 + t199) * t172 + (-t130 * t172 + (t131 - t211) * t173) * t146) * t165) * t137, 0, 0, -0.2e1 * t228 + 0.2e1 * (t130 * t146 * t137 + (-t137 * t227 - t146 * t228) * t148) * t148, 0;];
JaD_rot  = t1;
