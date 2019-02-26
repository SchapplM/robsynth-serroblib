% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRP3_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:02:23
% EndTime: 2019-02-26 20:02:23
% DurationCPUTime: 0.84s
% Computational Cost: add. (2688->101), mult. (8196->213), div. (535->12), fcn. (10572->13), ass. (0->99)
t187 = sin(pkin(10));
t190 = cos(pkin(10));
t193 = sin(qJ(2));
t191 = cos(pkin(6));
t195 = cos(qJ(2));
t218 = t191 * t195;
t176 = -t187 * t193 + t190 * t218;
t169 = t176 * qJD(2);
t219 = t191 * t193;
t177 = t187 * t195 + t190 * t219;
t192 = sin(qJ(3));
t188 = sin(pkin(6));
t222 = t188 * t192;
t209 = t190 * t222;
t194 = cos(qJ(3));
t215 = qJD(3) * t194;
t142 = -qJD(3) * t209 + t169 * t192 + t177 * t215;
t221 = t188 * t194;
t162 = t177 * t192 + t190 * t221;
t160 = t162 ^ 2;
t180 = -t191 * t194 + t193 * t222;
t174 = 0.1e1 / t180 ^ 2;
t156 = t160 * t174 + 0.1e1;
t154 = 0.1e1 / t156;
t181 = t191 * t192 + t193 * t221;
t216 = qJD(2) * t195;
t208 = t188 * t216;
t167 = qJD(3) * t181 + t192 * t208;
t173 = 0.1e1 / t180;
t226 = t162 * t174;
t126 = (-t142 * t173 + t167 * t226) * t154;
t157 = atan2(-t162, t180);
t152 = sin(t157);
t153 = cos(t157);
t206 = -t152 * t180 - t153 * t162;
t122 = t126 * t206 - t142 * t152 + t153 * t167;
t136 = -t152 * t162 + t153 * t180;
t133 = 0.1e1 / t136;
t134 = 0.1e1 / t136 ^ 2;
t239 = t122 * t133 * t134;
t210 = t187 * t219;
t179 = t190 * t195 - t210;
t166 = t179 * t194 + t187 * t222;
t178 = t187 * t218 + t190 * t193;
t186 = sin(pkin(11));
t189 = cos(pkin(11));
t150 = t166 * t186 - t178 * t189;
t171 = t178 * qJD(2);
t204 = -t179 * t192 + t187 * t221;
t145 = qJD(3) * t204 - t171 * t194;
t172 = -qJD(2) * t210 + t190 * t216;
t141 = t145 * t189 + t172 * t186;
t151 = t166 * t189 + t178 * t186;
t147 = 0.1e1 / t151;
t148 = 0.1e1 / t151 ^ 2;
t233 = t141 * t147 * t148;
t238 = 0.2e1 * t150 * t233;
t237 = -0.2e1 * t204 * t239;
t220 = t188 * t195;
t202 = -t173 * t176 + t220 * t226;
t236 = t192 * t202;
t225 = t167 * t173 * t174;
t235 = -0.2e1 * (t142 * t226 - t160 * t225) / t156 ^ 2;
t234 = t134 * t204;
t144 = qJD(3) * t166 - t171 * t192;
t232 = t144 * t134;
t231 = t147 * t186;
t230 = t148 * t150;
t229 = t150 * t189;
t228 = t152 * t204;
t227 = t153 * t204;
t224 = t178 * t192;
t223 = t178 * t194;
t217 = qJD(2) * t193;
t161 = t204 ^ 2;
t132 = t134 * t161 + 0.1e1;
t214 = 0.2e1 * (-t161 * t239 - t204 * t232) / t132 ^ 2;
t146 = t150 ^ 2;
t139 = t146 * t148 + 0.1e1;
t140 = t145 * t186 - t172 * t189;
t213 = 0.2e1 * (t140 * t230 - t146 * t233) / t139 ^ 2;
t207 = -0.2e1 * t162 * t225;
t164 = t177 * t194 - t209;
t205 = -t164 * t173 + t181 * t226;
t203 = qJD(3) * t224 - t172 * t194;
t170 = t177 * qJD(2);
t168 = -qJD(3) * t180 + t194 * t208;
t159 = t179 * t186 - t189 * t223;
t158 = -t179 * t189 - t186 * t223;
t143 = -qJD(3) * t162 + t169 * t194;
t137 = 0.1e1 / t139;
t129 = 0.1e1 / t132;
t128 = t154 * t236;
t127 = t205 * t154;
t124 = (-t152 * t176 + t153 * t220) * t192 + t206 * t128;
t123 = t127 * t206 - t152 * t164 + t153 * t181;
t121 = t205 * t235 + (t181 * t207 - t143 * t173 + (t142 * t181 + t162 * t168 + t164 * t167) * t174) * t154;
t119 = t235 * t236 + (t202 * t215 + (t207 * t220 + t170 * t173 + (t167 * t176 + (t142 * t195 - t162 * t217) * t188) * t174) * t192) * t154;
t1 = [0, t119, t121, 0, 0, 0; 0 (-t124 * t234 + t133 * t224) * t214 + ((-t172 * t192 - t178 * t215) * t133 + (-t232 + t237) * t124 + (t224 * t122 + (-t119 * t162 - t128 * t142 + (-t192 * t217 + t195 * t215) * t188 + (-t128 * t180 - t176 * t192) * t126) * t227 + (-t176 * t215 - t119 * t180 - t128 * t167 + t170 * t192 + (t128 * t162 - t192 * t220) * t126) * t228) * t134) * t129 (-t123 * t234 - t133 * t166) * t214 + (t123 * t237 + t145 * t133 + (-t166 * t122 - t123 * t144 + (-t121 * t162 - t127 * t142 + t168 + (-t127 * t180 - t164) * t126) * t227 + (-t121 * t180 - t127 * t167 - t143 + (t127 * t162 - t181) * t126) * t228) * t134) * t129, 0, 0, 0; 0 (-t147 * t158 + t159 * t230) * t213 + ((t171 * t189 + t186 * t203) * t147 + t159 * t238 + (-t158 * t141 - (-t171 * t186 + t189 * t203) * t150 - t159 * t140) * t148) * t137 -(-t229 * t148 + t231) * t204 * t213 + (t204 * t189 * t238 - t144 * t231 + (t144 * t229 - (t140 * t189 + t141 * t186) * t204) * t148) * t137, 0, 0, 0;];
JaD_rot  = t1;
