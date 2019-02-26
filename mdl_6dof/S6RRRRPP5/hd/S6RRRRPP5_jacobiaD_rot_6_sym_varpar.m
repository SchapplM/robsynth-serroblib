% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:27:39
% EndTime: 2019-02-26 22:27:40
% DurationCPUTime: 0.76s
% Computational Cost: add. (1236->97), mult. (2786->207), div. (506->12), fcn. (3255->9), ass. (0->99)
t157 = qJD(3) + qJD(4);
t165 = qJ(3) + qJ(4);
t155 = sin(t165);
t156 = cos(t165);
t167 = sin(qJ(1));
t168 = cos(qJ(2));
t169 = cos(qJ(1));
t205 = t168 * t169;
t231 = t155 * t205 - t167 * t156;
t236 = t231 * t157;
t166 = sin(qJ(2));
t159 = t166 ^ 2;
t162 = 0.1e1 / t168 ^ 2;
t213 = t159 * t162;
t160 = t167 ^ 2;
t153 = t160 * t213 + 0.1e1;
t161 = 0.1e1 / t168;
t210 = t161 * t166;
t233 = t166 * t213;
t177 = qJD(2) * (t161 * t233 + t210);
t202 = qJD(1) * t169;
t211 = t159 * t167;
t181 = t202 * t211;
t222 = (t160 * t177 + t162 * t181) / t153 ^ 2;
t235 = -0.2e1 * t222;
t199 = qJD(2) * t169;
t188 = t166 * t199;
t203 = qJD(1) * t168;
t125 = -t155 * t202 + (t167 * t203 + t188) * t156 + t236;
t139 = t231 ^ 2;
t146 = t167 * t155 + t156 * t205;
t141 = 0.1e1 / t146 ^ 2;
t220 = t139 * t141;
t132 = 0.1e1 + t220;
t183 = -t157 + t203;
t184 = t157 * t168 - qJD(1);
t214 = t156 * t169;
t124 = -t184 * t214 + (t183 * t167 + t188) * t155;
t217 = t141 * t231;
t196 = t124 * t217;
t140 = 0.1e1 / t146;
t142 = t140 * t141;
t219 = t139 * t142;
t234 = (t125 * t219 - t196) / t132 ^ 2;
t187 = 0.1e1 + t213;
t232 = t167 * t187;
t201 = qJD(2) * t167;
t230 = -t166 * t201 + t183 * t169;
t207 = t167 * t166;
t152 = atan2(t207, t168);
t148 = cos(t152);
t147 = sin(t152);
t193 = t147 * t207;
t136 = t148 * t168 + t193;
t133 = 0.1e1 / t136;
t134 = 0.1e1 / t136 ^ 2;
t229 = 0.2e1 * t166;
t150 = 0.1e1 / t153;
t228 = t150 - 0.1e1;
t164 = t169 ^ 2;
t212 = t159 * t164;
t131 = t134 * t212 + 0.1e1;
t200 = qJD(2) * t168;
t190 = t166 * t202;
t126 = ((t167 * t200 + t190) * t161 + t201 * t213) * t150;
t215 = t148 * t166;
t120 = (t126 * t167 - qJD(2)) * t215 + (t190 + (-t126 + t201) * t168) * t147;
t226 = t120 * t133 * t134;
t227 = (-t212 * t226 + (t164 * t166 * t200 - t181) * t134) / t131 ^ 2;
t225 = t126 * t147;
t224 = t126 * t166;
t223 = t134 * t166;
t138 = t150 * t232;
t221 = t138 * t167;
t218 = t140 * t155;
t216 = t231 * t156;
t209 = t166 * t169;
t206 = t167 * t168;
t204 = qJD(1) * t167;
t198 = 0.2e1 * t234;
t197 = -0.2e1 * t226;
t195 = t125 * t142 * t231;
t194 = t134 * t209;
t192 = t150 * t159 * t161;
t186 = -0.2e1 * t166 * t227;
t185 = t161 * t235;
t182 = t167 * t192;
t180 = t187 * t169;
t179 = t184 * t167;
t178 = -t141 * t216 + t218;
t144 = t155 * t169 - t156 * t206;
t143 = t155 * t206 + t214;
t129 = 0.1e1 / t132;
t127 = 0.1e1 / t131;
t123 = (-t228 * t166 * t147 + t148 * t182) * t169;
t122 = t147 * t206 - t215 + (-t147 * t168 + t148 * t207) * t138;
t121 = t232 * t235 + (qJD(1) * t180 + 0.2e1 * t167 * t177) * t150;
t117 = (t140 * t146 + t220) * t198 + (0.2e1 * t196 + (-t141 * t146 + t140 - 0.2e1 * t219) * t125) * t129;
t1 = [t185 * t209 + (qJD(2) * t180 - t204 * t210) * t150, t121, 0, 0, 0, 0; (t133 * t186 + (t133 * t200 + (-qJD(1) * t123 - t120) * t223) * t127) * t167 + (t134 * t186 * t123 + (((-t126 * t182 - t228 * t200 + t222 * t229) * t147 + (t185 * t211 + t224 + (-t224 + (t229 + t233) * t201) * t150) * t148) * t194 + (t134 * t200 + t166 * t197) * t123 + (t133 + ((-t160 + t164) * t148 * t192 + t228 * t193) * t134) * t166 * qJD(1)) * t127) * t169, 0.2e1 * (-t122 * t223 + t133 * t168) * t169 * t227 + ((t133 * t204 + (qJD(2) * t122 + t120) * t169 * t134) * t168 + (t133 * t199 + (t121 * t148 * t167 - t147 * t201 - t221 * t225 + t225 + (qJD(2) * t147 + t148 * t202) * t138) * t194 + (-t134 * t204 + t169 * t197) * t122 + ((-t121 + t202) * t147 + ((-0.1e1 + t221) * qJD(2) + (-t138 + t167) * t126) * t148) * t134 * t205) * t166) * t127, 0, 0, 0, 0; (-t140 * t143 - t144 * t217) * t198 + (0.2e1 * t144 * t195 + t140 * t156 * t179 + t230 * t218 + (t155 * t179 * t231 - t144 * t124 + t143 * t125 - t230 * t216) * t141) * t129, -0.2e1 * t178 * t209 * t234 + (t178 * t168 * t199 + (-t178 * t204 + ((t140 * t157 - 0.2e1 * t195) * t156 + (t124 * t156 + (t125 + t236) * t155) * t141) * t169) * t166) * t129, t117, t117, 0, 0;];
JaD_rot  = t1;
