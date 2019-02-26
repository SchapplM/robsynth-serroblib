% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR12_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:12
% EndTime: 2019-02-26 21:44:14
% DurationCPUTime: 0.66s
% Computational Cost: add. (1643->91), mult. (4303->199), div. (668->14), fcn. (5516->11), ass. (0->91)
t184 = sin(qJ(2));
t185 = sin(qJ(1));
t186 = cos(qJ(2));
t187 = cos(qJ(1));
t233 = cos(pkin(6));
t203 = t187 * t233;
t165 = t184 * t203 + t185 * t186;
t183 = sin(pkin(6));
t222 = t183 * t184;
t159 = atan2(-t165, t222);
t155 = sin(t159);
t156 = cos(t159);
t162 = t165 ^ 2;
t178 = 0.1e1 / t183 ^ 2;
t181 = 0.1e1 / t184 ^ 2;
t160 = t162 * t178 * t181 + 0.1e1;
t157 = 0.1e1 / t160;
t177 = 0.1e1 / t183;
t180 = 0.1e1 / t184;
t208 = t165 * t177 * t180;
t234 = (t156 * t208 + t155) * t157 - t155;
t139 = -t155 * t165 + t156 * t222;
t136 = 0.1e1 / t139;
t204 = t185 * t233;
t167 = t187 * t184 + t186 * t204;
t179 = qJ(4) + pkin(11);
t175 = sin(t179);
t176 = cos(t179);
t221 = t183 * t185;
t152 = t167 * t175 + t176 * t221;
t148 = 0.1e1 / t152;
t137 = 0.1e1 / t139 ^ 2;
t149 = 0.1e1 / t152 ^ 2;
t199 = qJD(2) * t233 + qJD(1);
t201 = t184 * t204;
t216 = qJD(2) * t184;
t218 = t187 * t186;
t146 = -qJD(1) * t201 - t185 * t216 + t199 * t218;
t215 = qJD(2) * t186;
t205 = t181 * t215;
t194 = -t146 * t180 + t165 * t205;
t224 = t157 * t177;
t128 = t194 * t224;
t196 = -t155 * t222 - t156 * t165;
t209 = t156 * t183 * t186;
t124 = qJD(2) * t209 + t196 * t128 - t155 * t146;
t232 = t124 * t136 * t137;
t200 = t186 * t203;
t219 = t185 * t184;
t143 = -qJD(1) * t200 - t187 * t215 + t199 * t219;
t217 = qJD(1) * t183;
t206 = t187 * t217;
t151 = -t167 * t176 + t175 * t221;
t214 = qJD(4) * t151;
t134 = -t143 * t175 + t176 * t206 - t214;
t231 = t134 * t148 * t149;
t182 = t180 * t181;
t230 = (t146 * t165 * t181 - t162 * t182 * t215) * t178 / t160 ^ 2;
t168 = -t201 + t218;
t229 = t137 * t168;
t133 = t152 * qJD(4) + t143 * t176 + t175 * t206;
t147 = t151 ^ 2;
t142 = t147 * t149 + 0.1e1;
t227 = t149 * t151;
t228 = 0.1e1 / t142 ^ 2 * (t133 * t227 - t147 * t231);
t226 = t155 * t168;
t225 = t156 * t168;
t223 = t181 * t186;
t220 = t183 * t187;
t163 = t168 ^ 2;
t132 = t137 * t163 + 0.1e1;
t144 = t165 * qJD(1) + t167 * qJD(2);
t213 = 0.2e1 * (-t144 * t229 - t163 * t232) / t132 ^ 2;
t212 = 0.2e1 * t232;
t211 = -0.2e1 * t230;
t210 = 0.2e1 * t228;
t207 = t185 * t217;
t202 = 0.2e1 * t151 * t231;
t197 = t148 * t176 + t175 * t227;
t164 = -t200 + t219;
t195 = t164 * t180 + t165 * t223;
t154 = -t164 * t175 + t176 * t220;
t153 = t164 * t176 + t175 * t220;
t145 = t167 * qJD(1) + t165 * qJD(2);
t140 = 0.1e1 / t142;
t130 = 0.1e1 / t132;
t129 = t195 * t224;
t127 = t234 * t168;
t125 = t196 * t129 + t155 * t164 + t209;
t123 = (t195 * t211 + (t146 * t223 + t145 * t180 + (-t164 * t223 + (-0.2e1 * t182 * t186 ^ 2 - t180) * t165) * qJD(2)) * t157) * t177;
t1 = [(0.2e1 * t168 * t180 * t230 + (t144 * t180 + t168 * t205) * t157) * t177, t123, 0, 0, 0, 0; t165 * t136 * t213 + (-t146 * t136 + (t124 * t165 + t127 * t144) * t137) * t130 + (t127 * t212 * t130 + (t127 * t213 + (-(-t128 * t157 * t208 + t211) * t226 - (t208 * t211 - t128 + (-t194 * t177 + t128) * t157) * t225 + t234 * t144) * t130) * t137) * t168 (t125 * t229 + t136 * t167) * t213 + (t125 * t168 * t212 + t143 * t136 + (t167 * t124 + t125 * t144 - (-t183 * t216 - t123 * t165 - t129 * t146 + (-t129 * t222 + t164) * t128) * t225 - (t128 * t129 * t165 + t145 + (-t123 * t184 + (-qJD(2) * t129 - t128) * t186) * t183) * t226) * t137) * t130, 0, 0, 0, 0; (-t148 * t153 + t154 * t227) * t210 + ((t154 * qJD(4) + t145 * t176 - t175 * t207) * t148 + t154 * t202 + (-t153 * t134 - (-t153 * qJD(4) - t145 * t175 - t176 * t207) * t151 - t154 * t133) * t149) * t140, t197 * t168 * t210 + (t197 * t144 + ((qJD(4) * t148 + t202) * t175 + (-t133 * t175 + (t134 - t214) * t176) * t149) * t168) * t140, 0, -0.2e1 * t228 + 0.2e1 * (t133 * t149 * t140 + (-t140 * t231 - t149 * t228) * t151) * t151, 0, 0;];
JaD_rot  = t1;
