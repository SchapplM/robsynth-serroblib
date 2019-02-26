% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:10
% EndTime: 2019-02-26 21:55:11
% DurationCPUTime: 0.72s
% Computational Cost: add. (2887->96), mult. (2734->202), div. (498->12), fcn. (3199->9), ass. (0->96)
t165 = qJ(2) + pkin(11);
t160 = sin(t165);
t156 = t160 ^ 2;
t161 = cos(t165);
t158 = 0.1e1 / t161 ^ 2;
t214 = t156 * t158;
t169 = sin(qJ(1));
t232 = 0.2e1 * t169;
t231 = t160 * t214;
t166 = t169 ^ 2;
t152 = t166 * t214 + 0.1e1;
t150 = 0.1e1 / t152;
t157 = 0.1e1 / t161;
t170 = cos(qJ(1));
t204 = qJD(1) * t170;
t192 = t160 * t204;
t202 = qJD(2) * t169;
t123 = (-(-t161 * t202 - t192) * t157 + t202 * t214) * t150;
t230 = t123 - t202;
t168 = qJ(4) + qJ(5);
t163 = cos(t168);
t206 = t170 * t163;
t162 = sin(t168);
t209 = t169 * t162;
t145 = t161 * t206 + t209;
t210 = t169 * t160;
t148 = atan2(-t210, -t161);
t147 = cos(t148);
t146 = sin(t148);
t195 = t146 * t210;
t132 = -t147 * t161 - t195;
t129 = 0.1e1 / t132;
t139 = 0.1e1 / t145;
t130 = 0.1e1 / t132 ^ 2;
t140 = 0.1e1 / t145 ^ 2;
t229 = -0.2e1 * t160;
t228 = t150 - 0.1e1;
t216 = t147 * t160;
t118 = (-t123 * t169 + qJD(2)) * t216 + (t230 * t161 - t192) * t146;
t227 = t118 * t129 * t130;
t164 = qJD(4) + qJD(5);
t180 = t161 * t209 + t206;
t201 = qJD(2) * t170;
t191 = t160 * t201;
t124 = t180 * qJD(1) - t145 * t164 + t162 * t191;
t207 = t170 * t162;
t208 = t169 * t163;
t144 = t161 * t207 - t208;
t138 = t144 ^ 2;
t137 = t138 * t140 + 0.1e1;
t219 = t140 * t144;
t185 = -qJD(1) * t161 + t164;
t186 = t161 * t164 - qJD(1);
t125 = -t186 * t207 + (t185 * t169 - t191) * t163;
t224 = t125 * t139 * t140;
t226 = (-t124 * t219 - t138 * t224) / t137 ^ 2;
t225 = t123 * t160;
t223 = t130 * t160;
t222 = t130 * t170;
t212 = t157 * t160;
t179 = qJD(2) * (t157 * t231 + t212);
t183 = t156 * t169 * t204;
t221 = (t158 * t183 + t166 * t179) / t152 ^ 2;
t220 = t139 * t162;
t218 = t144 * t163;
t217 = t146 * t169;
t215 = t156 * t157;
t167 = t170 ^ 2;
t213 = t156 * t167;
t211 = t160 * t170;
t205 = qJD(1) * t169;
t203 = qJD(2) * t161;
t128 = t130 * t213 + 0.1e1;
t200 = 0.2e1 * (-t213 * t227 + (t160 * t167 * t203 - t183) * t130) / t128 ^ 2;
t199 = 0.2e1 * t227;
t198 = -0.2e1 * t226;
t197 = t130 * t211;
t196 = t144 * t224;
t194 = t150 * t215;
t190 = 0.1e1 + t214;
t189 = t160 * t200;
t188 = t221 * t229;
t187 = t221 * t232;
t184 = t169 * t194;
t182 = t190 * t170;
t181 = t140 * t218 - t220;
t178 = t160 * t202 + t185 * t170;
t143 = -t161 * t208 + t207;
t135 = 0.1e1 / t137;
t134 = t190 * t169 * t150;
t126 = 0.1e1 / t128;
t122 = (t228 * t160 * t146 - t147 * t184) * t170;
t121 = -t161 * t217 + t216 + (t146 * t161 - t147 * t210) * t134;
t119 = -t190 * t187 + (qJD(1) * t182 + t179 * t232) * t150;
t116 = t198 + 0.2e1 * (-t124 * t140 * t135 + (-t135 * t224 - t140 * t226) * t144) * t144;
t1 = [t170 * t157 * t188 + (qJD(2) * t182 - t205 * t212) * t150, t119, 0, 0, 0, 0; (t129 * t189 + (-t129 * t203 + (qJD(1) * t122 + t118) * t223) * t126) * t169 + (t130 * t189 * t122 + (-((t123 * t184 + t228 * t203 + t188) * t146 + (t187 * t215 - t225 + (t225 + (t229 - t231) * t202) * t150) * t147) * t197 + (-t130 * t203 + t160 * t199) * t122 + (-t129 + ((-t166 + t167) * t147 * t194 + t228 * t195) * t130) * t160 * qJD(1)) * t126) * t170 (t121 * t223 - t129 * t161) * t170 * t200 + ((-t129 * t205 + (-qJD(2) * t121 - t118) * t222) * t161 + (-t129 * t201 - (-t119 * t147 * t169 - t230 * t146 + (-qJD(2) * t146 + t123 * t217 - t147 * t204) * t134) * t197 + (t130 * t205 + t170 * t199) * t121 - ((t119 - t204) * t146 + ((-t134 * t169 + 0.1e1) * qJD(2) + (t134 - t169) * t123) * t147) * t161 * t222) * t160) * t126, 0, 0, 0, 0; 0.2e1 * (t139 * t180 + t143 * t219) * t226 + (0.2e1 * t143 * t196 - t186 * t139 * t208 + t178 * t220 + (-t144 * t186 * t209 + t143 * t124 + t125 * t180 - t178 * t218) * t140) * t135, t181 * t198 * t211 + (t181 * t161 * t201 + (-t181 * t205 + ((-t139 * t164 - 0.2e1 * t196) * t163 + (-t124 * t163 + (-t144 * t164 + t125) * t162) * t140) * t170) * t160) * t135, 0, t116, t116, 0;];
JaD_rot  = t1;
