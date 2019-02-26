% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:28:45
% EndTime: 2019-02-26 21:28:45
% DurationCPUTime: 0.72s
% Computational Cost: add. (3326->96), mult. (2734->202), div. (498->12), fcn. (3199->9), ass. (0->96)
t166 = qJ(2) + pkin(10);
t162 = sin(t166);
t158 = t162 ^ 2;
t163 = cos(t166);
t160 = 0.1e1 / t163 ^ 2;
t214 = t158 * t160;
t169 = sin(qJ(1));
t232 = 0.2e1 * t169;
t231 = t162 * t214;
t167 = t169 ^ 2;
t152 = t167 * t214 + 0.1e1;
t150 = 0.1e1 / t152;
t159 = 0.1e1 / t163;
t170 = cos(qJ(1));
t204 = qJD(1) * t170;
t192 = t162 * t204;
t202 = qJD(2) * t169;
t125 = (-(-t163 * t202 - t192) * t159 + t202 * t214) * t150;
t230 = t125 - t202;
t164 = pkin(11) + qJ(5) + qJ(6);
t156 = cos(t164);
t206 = t170 * t156;
t155 = sin(t164);
t210 = t169 * t155;
t145 = t163 * t206 + t210;
t208 = t169 * t162;
t149 = atan2(-t208, -t163);
t148 = cos(t149);
t147 = sin(t149);
t195 = t147 * t208;
t135 = -t148 * t163 - t195;
t132 = 0.1e1 / t135;
t139 = 0.1e1 / t145;
t133 = 0.1e1 / t135 ^ 2;
t140 = 0.1e1 / t145 ^ 2;
t229 = -0.2e1 * t162;
t228 = t150 - 0.1e1;
t216 = t148 * t162;
t118 = (-t125 * t169 + qJD(2)) * t216 + (t230 * t163 - t192) * t147;
t227 = t118 * t132 * t133;
t165 = qJD(5) + qJD(6);
t180 = t163 * t210 + t206;
t201 = qJD(2) * t170;
t191 = t162 * t201;
t123 = qJD(1) * t180 - t145 * t165 + t155 * t191;
t207 = t170 * t155;
t209 = t169 * t156;
t144 = t163 * t207 - t209;
t138 = t144 ^ 2;
t131 = t138 * t140 + 0.1e1;
t219 = t140 * t144;
t185 = -qJD(1) * t163 + t165;
t186 = t163 * t165 - qJD(1);
t124 = -t186 * t207 + (t169 * t185 - t191) * t156;
t225 = t124 * t139 * t140;
t226 = (-t123 * t219 - t138 * t225) / t131 ^ 2;
t224 = t125 * t162;
t223 = t133 * t162;
t222 = t133 * t170;
t212 = t159 * t162;
t179 = qJD(2) * (t159 * t231 + t212);
t183 = t158 * t169 * t204;
t221 = (t160 * t183 + t167 * t179) / t152 ^ 2;
t220 = t139 * t155;
t218 = t144 * t156;
t217 = t147 * t169;
t215 = t158 * t159;
t168 = t170 ^ 2;
t213 = t158 * t168;
t211 = t162 * t170;
t205 = qJD(1) * t169;
t203 = qJD(2) * t163;
t128 = t133 * t213 + 0.1e1;
t200 = 0.2e1 * (-t213 * t227 + (t162 * t168 * t203 - t183) * t133) / t128 ^ 2;
t199 = 0.2e1 * t227;
t198 = -0.2e1 * t226;
t197 = t133 * t211;
t196 = t144 * t225;
t194 = t150 * t215;
t190 = 0.1e1 + t214;
t189 = t162 * t200;
t188 = t221 * t229;
t187 = t221 * t232;
t184 = t169 * t194;
t182 = t190 * t170;
t181 = t140 * t218 - t220;
t178 = t162 * t202 + t170 * t185;
t143 = -t163 * t209 + t207;
t137 = t190 * t169 * t150;
t129 = 0.1e1 / t131;
t126 = 0.1e1 / t128;
t122 = (t228 * t162 * t147 - t148 * t184) * t170;
t121 = -t163 * t217 + t216 + (t147 * t163 - t148 * t208) * t137;
t119 = -t190 * t187 + (qJD(1) * t182 + t179 * t232) * t150;
t116 = t198 + 0.2e1 * (-t123 * t140 * t129 + (-t129 * t225 - t140 * t226) * t144) * t144;
t1 = [t170 * t159 * t188 + (qJD(2) * t182 - t205 * t212) * t150, t119, 0, 0, 0, 0; (t132 * t189 + (-t132 * t203 + (qJD(1) * t122 + t118) * t223) * t126) * t169 + (t133 * t189 * t122 + (-((t125 * t184 + t228 * t203 + t188) * t147 + (t187 * t215 - t224 + (t224 + (t229 - t231) * t202) * t150) * t148) * t197 + (-t133 * t203 + t162 * t199) * t122 + (-t132 + ((-t167 + t168) * t148 * t194 + t228 * t195) * t133) * t162 * qJD(1)) * t126) * t170 (t121 * t223 - t132 * t163) * t170 * t200 + ((-t132 * t205 + (-qJD(2) * t121 - t118) * t222) * t163 + (-t132 * t201 - (-t119 * t148 * t169 - t230 * t147 + (-qJD(2) * t147 + t125 * t217 - t148 * t204) * t137) * t197 + (t133 * t205 + t170 * t199) * t121 - ((t119 - t204) * t147 + ((-t137 * t169 + 0.1e1) * qJD(2) + (t137 - t169) * t125) * t148) * t163 * t222) * t162) * t126, 0, 0, 0, 0; 0.2e1 * (t139 * t180 + t143 * t219) * t226 + (0.2e1 * t143 * t196 - t186 * t139 * t209 + t178 * t220 + (-t144 * t186 * t210 + t143 * t123 + t124 * t180 - t178 * t218) * t140) * t129, t181 * t198 * t211 + (t181 * t163 * t201 + (-t181 * t205 + ((-t139 * t165 - 0.2e1 * t196) * t156 + (-t123 * t156 + (-t144 * t165 + t124) * t155) * t140) * t170) * t162) * t129, 0, 0, t116, t116;];
JaD_rot  = t1;
