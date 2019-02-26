% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP8_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:11:54
% EndTime: 2019-02-26 21:11:55
% DurationCPUTime: 0.71s
% Computational Cost: add. (3360->95), mult. (3810->207), div. (753->12), fcn. (4455->9), ass. (0->97)
t168 = cos(qJ(1));
t228 = 0.2e1 * t168;
t163 = t168 ^ 2;
t164 = qJ(3) + qJ(4);
t159 = sin(t164);
t155 = 0.1e1 / t159 ^ 2;
t160 = cos(t164);
t158 = t160 ^ 2;
t211 = t155 * t158;
t152 = t163 * t211 + 0.1e1;
t150 = 0.1e1 / t152;
t154 = 0.1e1 / t159;
t166 = sin(qJ(1));
t200 = qJD(1) * t166;
t192 = t160 * t200;
t161 = qJD(3) + qJD(4);
t206 = t161 * t168;
t193 = t155 * t206;
t124 = ((t159 * t206 + t192) * t154 + t158 * t193) * t150;
t227 = -t124 + t206;
t202 = t168 * t160;
t149 = atan2(-t202, t159);
t147 = sin(t149);
t148 = cos(t149);
t134 = -t147 * t202 + t148 * t159;
t131 = 0.1e1 / t134;
t167 = cos(qJ(5));
t203 = t166 * t167;
t165 = sin(qJ(5));
t205 = t165 * t168;
t144 = t159 * t203 + t205;
t140 = 0.1e1 / t144;
t132 = 0.1e1 / t134 ^ 2;
t141 = 0.1e1 / t144 ^ 2;
t162 = t166 ^ 2;
t210 = t158 * t162;
t127 = t132 * t210 + 0.1e1;
t199 = qJD(1) * t168;
t184 = t158 * t166 * t199;
t208 = t160 * t161;
t214 = t148 * t160;
t223 = t124 * t168;
t119 = (t161 - t223) * t214 + (t159 * t227 + t192) * t147;
t225 = t119 * t131 * t132;
t226 = (-t210 * t225 + (-t159 * t162 * t208 + t184) * t132) / t127 ^ 2;
t187 = qJD(1) * t159 + qJD(5);
t180 = t187 * t168;
t188 = qJD(5) * t159 + qJD(1);
t181 = t188 * t167;
t207 = t161 * t166;
t128 = t166 * t181 + (t160 * t207 + t180) * t165;
t201 = t168 * t167;
t204 = t166 * t165;
t143 = t159 * t204 - t201;
t139 = t143 ^ 2;
t138 = t139 * t141 + 0.1e1;
t217 = t141 * t143;
t182 = t188 * t165;
t129 = t167 * t180 + (t167 * t208 - t182) * t166;
t222 = t129 * t140 * t141;
t224 = (t128 * t217 - t139 * t222) / t138 ^ 2;
t157 = t160 * t158;
t212 = t154 * t160;
t178 = t161 * (-t154 * t155 * t157 - t212);
t221 = (-t155 * t184 + t163 * t178) / t152 ^ 2;
t220 = t132 * t160;
t219 = t132 * t166;
t218 = t140 * t165;
t216 = t143 * t167;
t215 = t147 * t168;
t213 = t154 * t158;
t209 = t159 * t161;
t198 = -0.2e1 * t225;
t197 = 0.2e1 * t224;
t196 = t160 * t226;
t195 = t160 * t221;
t194 = t160 * t219;
t191 = 0.1e1 + t211;
t190 = 0.2e1 * t143 * t222;
t189 = t221 * t228;
t186 = t148 * t150 * t213;
t185 = (-t150 + 0.1e1) * t160 * t147;
t183 = t191 * t166;
t179 = t141 * t216 - t218;
t177 = t179 * t166;
t176 = t161 * t202 - t187 * t166;
t146 = t159 * t201 - t204;
t145 = t159 * t205 + t203;
t136 = 0.1e1 / t138;
t135 = t191 * t168 * t150;
t125 = 0.1e1 / t127;
t123 = (-t168 * t186 + t185) * t166;
t121 = t159 * t215 + t214 + (-t147 * t159 - t148 * t202) * t135;
t120 = -t191 * t189 + (-qJD(1) * t183 + t178 * t228) * t150;
t117 = t160 * t177 * t197 + (t177 * t209 + (-t179 * t199 + ((qJD(5) * t140 + t190) * t167 + (-t128 * t167 + (qJD(5) * t143 - t129) * t165) * t141) * t166) * t160) * t136;
t116 = 0.2e1 * (-t121 * t220 - t131 * t159) * t166 * t226 + ((t131 * t199 + (-t121 * t161 - t119) * t219) * t159 + (t131 * t207 + (-t120 * t148 * t168 + t227 * t147 + (t124 * t215 - t147 * t161 + t148 * t200) * t135) * t194 + (t132 * t199 + t166 * t198) * t121 + ((-t120 - t200) * t147 + ((t135 * t168 - 0.1e1) * t161 + (-t135 + t168) * t124) * t148) * t159 * t219) * t160) * t125;
t1 = [-0.2e1 * t166 * t154 * t195 + (-t161 * t183 + t199 * t212) * t150, 0, t120, t120, 0, 0; (0.2e1 * t131 * t196 + (t131 * t209 + (qJD(1) * t123 + t119) * t220) * t125) * t168 + (-0.2e1 * t132 * t196 * t123 + (((0.2e1 * t195 - t209 + (t213 * t223 + t209) * t150) * t147 + (t189 * t213 + t124 * t160 + (t157 * t193 + (-t124 + 0.2e1 * t206) * t160) * t150) * t148) * t194 + (-t132 * t209 + t160 * t198) * t123 + (t131 + ((t162 - t163) * t186 + t168 * t185) * t132) * t160 * qJD(1)) * t125) * t166, 0, t116, t116, 0, 0; (-t140 * t145 + t146 * t217) * t197 + (t146 * t190 + t168 * t140 * t181 + t176 * t218 + (t168 * t143 * t182 - t146 * t128 - t145 * t129 - t176 * t216) * t141) * t136, 0, t117, t117, -0.2e1 * t224 + 0.2e1 * (t128 * t136 * t141 + (-t136 * t222 - t141 * t224) * t143) * t143, 0;];
JaD_rot  = t1;
