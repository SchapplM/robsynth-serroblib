% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPPR5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:00:42
% EndTime: 2019-02-26 20:00:43
% DurationCPUTime: 0.84s
% Computational Cost: add. (2688->99), mult. (8196->209), div. (535->12), fcn. (10572->13), ass. (0->99)
t190 = cos(pkin(6));
t194 = cos(qJ(2));
t239 = cos(pkin(10));
t212 = t239 * t194;
t187 = sin(pkin(10));
t192 = sin(qJ(2));
t226 = t187 * t192;
t177 = t190 * t212 - t226;
t170 = t177 * qJD(2);
t193 = cos(qJ(3));
t191 = sin(qJ(3));
t213 = t239 * t192;
t225 = t187 * t194;
t202 = -t190 * t213 - t225;
t188 = sin(pkin(6));
t214 = t188 * t239;
t241 = t202 * t191 - t193 * t214;
t143 = t241 * qJD(3) + t170 * t193;
t163 = -t191 * t214 - t202 * t193;
t160 = t163 ^ 2;
t223 = t188 * t193;
t181 = t190 * t191 + t192 * t223;
t175 = 0.1e1 / t181 ^ 2;
t156 = t160 * t175 + 0.1e1;
t154 = 0.1e1 / t156;
t224 = t188 * t191;
t180 = t190 * t193 - t192 * t224;
t222 = t188 * t194;
t215 = qJD(2) * t222;
t168 = t180 * qJD(3) + t193 * t215;
t174 = 0.1e1 / t181;
t230 = t163 * t175;
t126 = (-t143 * t174 + t168 * t230) * t154;
t157 = atan2(-t163, t181);
t152 = sin(t157);
t153 = cos(t157);
t208 = -t152 * t181 - t153 * t163;
t122 = t208 * t126 - t143 * t152 + t153 * t168;
t136 = -t152 * t163 + t153 * t181;
t133 = 0.1e1 / t136;
t134 = 0.1e1 / t136 ^ 2;
t244 = t122 * t133 * t134;
t179 = -t190 * t226 + t212;
t166 = t179 * t193 + t187 * t224;
t243 = 0.2e1 * t166 * t244;
t204 = -t174 * t177 + t222 * t230;
t242 = t193 * t204;
t229 = t168 * t174 * t175;
t240 = -0.2e1 * (t143 * t230 - t160 * t229) / t156 ^ 2;
t186 = sin(pkin(11));
t189 = cos(pkin(11));
t203 = -t190 * t225 - t213;
t206 = -t179 * t191 + t187 * t223;
t151 = -t186 * t206 - t189 * t203;
t147 = 0.1e1 / t151;
t148 = 0.1e1 / t151 ^ 2;
t238 = t134 * t166;
t172 = t203 * qJD(2);
t144 = t166 * qJD(3) + t172 * t191;
t173 = t179 * qJD(2);
t141 = t144 * t186 + t173 * t189;
t237 = t141 * t147 * t148;
t145 = t206 * qJD(3) + t172 * t193;
t236 = t145 * t134;
t235 = t147 * t189;
t150 = -t186 * t203 + t189 * t206;
t234 = t148 * t150;
t233 = t150 * t186;
t232 = t152 * t166;
t231 = t153 * t166;
t228 = t203 * t191;
t227 = t203 * t193;
t221 = qJD(2) * t192;
t220 = qJD(3) * t191;
t161 = t166 ^ 2;
t132 = t134 * t161 + 0.1e1;
t219 = 0.2e1 * (-t161 * t244 + t166 * t236) / t132 ^ 2;
t146 = t150 ^ 2;
t139 = t146 * t148 + 0.1e1;
t140 = -t144 * t189 + t173 * t186;
t218 = 0.2e1 * (t140 * t234 - t146 * t237) / t139 ^ 2;
t211 = 0.2e1 * t150 * t237;
t210 = -0.2e1 * t163 * t229;
t207 = -t174 * t241 + t180 * t230;
t205 = -qJD(3) * t227 + t173 * t191;
t171 = t202 * qJD(2);
t167 = -t181 * qJD(3) - t191 * t215;
t159 = t179 * t189 + t186 * t228;
t158 = t179 * t186 - t189 * t228;
t142 = t163 * qJD(3) + t170 * t191;
t137 = 0.1e1 / t139;
t129 = 0.1e1 / t132;
t128 = t154 * t242;
t127 = t207 * t154;
t124 = (-t152 * t177 + t153 * t222) * t193 + t208 * t128;
t123 = t208 * t127 - t152 * t241 + t153 * t180;
t121 = t207 * t240 + (t180 * t210 + t142 * t174 + (t143 * t180 + t163 * t167 + t168 * t241) * t175) * t154;
t119 = t240 * t242 + (-t204 * t220 + (t210 * t222 - t171 * t174 + (t168 * t177 + (t143 * t194 - t163 * t221) * t188) * t175) * t193) * t154;
t1 = [0, t119, t121, 0, 0, 0; 0 (t124 * t238 - t133 * t227) * t219 + ((-t173 * t193 - t203 * t220) * t133 + (-t236 + t243) * t124 + (-t227 * t122 - (-t119 * t163 - t128 * t143 + (-t193 * t221 - t194 * t220) * t188 + (-t128 * t181 - t177 * t193) * t126) * t231 - (t177 * t220 - t119 * t181 - t128 * t168 - t171 * t193 + (t128 * t163 - t193 * t222) * t126) * t232) * t134) * t129 (t123 * t238 - t133 * t206) * t219 + (t123 * t243 - t144 * t133 + (-t206 * t122 - t123 * t145 - (-t121 * t163 - t127 * t143 + t167 + (-t127 * t181 - t241) * t126) * t231 - (-t121 * t181 - t127 * t168 + t142 + (t127 * t163 - t180) * t126) * t232) * t134) * t129, 0, 0, 0; 0 (-t147 * t158 + t159 * t234) * t218 + ((t172 * t186 + t205 * t189) * t147 + t159 * t211 + (-t158 * t141 - (t172 * t189 - t205 * t186) * t150 - t159 * t140) * t148) * t137 (t148 * t233 + t235) * t166 * t218 + (t166 * t186 * t211 - t145 * t235 + (-t145 * t233 + (-t140 * t186 + t141 * t189) * t166) * t148) * t137, 0, 0, 0;];
JaD_rot  = t1;
