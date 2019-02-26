% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:10:50
% EndTime: 2019-02-26 21:10:51
% DurationCPUTime: 0.79s
% Computational Cost: add. (2887->96), mult. (2734->202), div. (498->12), fcn. (3199->9), ass. (0->96)
t162 = pkin(10) + qJ(3);
t158 = sin(t162);
t154 = t158 ^ 2;
t159 = cos(t162);
t156 = 0.1e1 / t159 ^ 2;
t212 = t154 * t156;
t167 = sin(qJ(1));
t230 = 0.2e1 * t167;
t229 = t158 * t212;
t164 = t167 ^ 2;
t150 = t164 * t212 + 0.1e1;
t148 = 0.1e1 / t150;
t155 = 0.1e1 / t159;
t168 = cos(qJ(1));
t202 = qJD(1) * t168;
t190 = t158 * t202;
t200 = qJD(3) * t167;
t121 = (-(-t200 * t159 - t190) * t155 + t200 * t212) * t148;
t228 = t121 - t200;
t166 = qJ(4) + qJ(5);
t161 = cos(t166);
t204 = t168 * t161;
t160 = sin(t166);
t207 = t167 * t160;
t143 = t159 * t204 + t207;
t208 = t167 * t158;
t146 = atan2(-t208, -t159);
t145 = cos(t146);
t144 = sin(t146);
t193 = t144 * t208;
t130 = -t145 * t159 - t193;
t127 = 0.1e1 / t130;
t137 = 0.1e1 / t143;
t128 = 0.1e1 / t130 ^ 2;
t138 = 0.1e1 / t143 ^ 2;
t227 = -0.2e1 * t158;
t226 = t148 - 0.1e1;
t214 = t145 * t158;
t116 = (-t121 * t167 + qJD(3)) * t214 + (t228 * t159 - t190) * t144;
t225 = t116 * t127 * t128;
t163 = qJD(4) + qJD(5);
t178 = t159 * t207 + t204;
t199 = qJD(3) * t168;
t189 = t158 * t199;
t122 = qJD(1) * t178 - t143 * t163 + t160 * t189;
t205 = t168 * t160;
t206 = t167 * t161;
t142 = t159 * t205 - t206;
t136 = t142 ^ 2;
t135 = t136 * t138 + 0.1e1;
t217 = t138 * t142;
t183 = -qJD(1) * t159 + t163;
t184 = t159 * t163 - qJD(1);
t123 = -t184 * t205 + (t167 * t183 - t189) * t161;
t222 = t123 * t137 * t138;
t224 = (-t122 * t217 - t136 * t222) / t135 ^ 2;
t223 = t121 * t158;
t221 = t128 * t158;
t220 = t128 * t168;
t210 = t155 * t158;
t177 = qJD(3) * (t155 * t229 + t210);
t181 = t154 * t167 * t202;
t219 = (t156 * t181 + t164 * t177) / t150 ^ 2;
t218 = t137 * t160;
t216 = t142 * t161;
t215 = t144 * t167;
t213 = t154 * t155;
t165 = t168 ^ 2;
t211 = t154 * t165;
t209 = t158 * t168;
t203 = qJD(1) * t167;
t201 = qJD(3) * t159;
t126 = t128 * t211 + 0.1e1;
t198 = 0.2e1 * (-t211 * t225 + (t158 * t165 * t201 - t181) * t128) / t126 ^ 2;
t197 = 0.2e1 * t225;
t196 = -0.2e1 * t224;
t195 = t128 * t209;
t194 = t142 * t222;
t192 = t148 * t213;
t188 = 0.1e1 + t212;
t187 = t158 * t198;
t186 = t219 * t227;
t185 = t219 * t230;
t182 = t167 * t192;
t180 = t188 * t168;
t179 = t216 * t138 - t218;
t176 = t158 * t200 + t168 * t183;
t141 = -t159 * t206 + t205;
t133 = 0.1e1 / t135;
t132 = t188 * t167 * t148;
t124 = 0.1e1 / t126;
t120 = (t144 * t158 * t226 - t145 * t182) * t168;
t119 = -t159 * t215 + t214 + (t144 * t159 - t145 * t208) * t132;
t117 = -t188 * t185 + (qJD(1) * t180 + t177 * t230) * t148;
t114 = t196 + 0.2e1 * (-t122 * t138 * t133 + (-t133 * t222 - t138 * t224) * t142) * t142;
t1 = [t168 * t155 * t186 + (qJD(3) * t180 - t203 * t210) * t148, 0, t117, 0, 0, 0; (t127 * t187 + (-t127 * t201 + (qJD(1) * t120 + t116) * t221) * t124) * t167 + (t128 * t187 * t120 + (-((t121 * t182 + t201 * t226 + t186) * t144 + (t185 * t213 - t223 + (t223 + (t227 - t229) * t200) * t148) * t145) * t195 + (-t128 * t201 + t158 * t197) * t120 + (-t127 + ((-t164 + t165) * t145 * t192 + t226 * t193) * t128) * t158 * qJD(1)) * t124) * t168, 0 (t119 * t221 - t127 * t159) * t168 * t198 + ((-t127 * t203 + (-qJD(3) * t119 - t116) * t220) * t159 + (-t127 * t199 - (-t117 * t145 * t167 - t228 * t144 + (-qJD(3) * t144 + t121 * t215 - t145 * t202) * t132) * t195 + (t128 * t203 + t168 * t197) * t119 - ((t117 - t202) * t144 + ((-t132 * t167 + 0.1e1) * qJD(3) + (t132 - t167) * t121) * t145) * t159 * t220) * t158) * t124, 0, 0, 0; 0.2e1 * (t137 * t178 + t141 * t217) * t224 + (0.2e1 * t141 * t194 - t184 * t137 * t206 + t176 * t218 + (-t142 * t184 * t207 + t141 * t122 + t123 * t178 - t176 * t216) * t138) * t133, 0, t179 * t196 * t209 + (t179 * t159 * t199 + (-t179 * t203 + ((-t137 * t163 - 0.2e1 * t194) * t161 + (-t122 * t161 + (-t142 * t163 + t123) * t160) * t138) * t168) * t158) * t133, t114, t114, 0;];
JaD_rot  = t1;
