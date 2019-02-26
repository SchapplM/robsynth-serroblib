% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:52:47
% EndTime: 2019-02-26 20:52:47
% DurationCPUTime: 0.71s
% Computational Cost: add. (2698->94), mult. (2734->203), div. (498->12), fcn. (3199->9), ass. (0->97)
t162 = qJ(3) + pkin(10);
t157 = sin(t162);
t153 = 0.1e1 / t157 ^ 2;
t158 = cos(t162);
t156 = t158 ^ 2;
t209 = t153 * t156;
t167 = cos(qJ(1));
t228 = 0.2e1 * t167;
t227 = t158 * t209;
t164 = t167 ^ 2;
t150 = t164 * t209 + 0.1e1;
t148 = 0.1e1 / t150;
t152 = 0.1e1 / t157;
t166 = sin(qJ(1));
t202 = qJD(1) * t166;
t191 = t158 * t202;
t198 = qJD(3) * t167;
t122 = ((t157 * t198 + t191) * t152 + t198 * t209) * t148;
t226 = -t122 + t198;
t204 = t167 * t158;
t147 = atan2(-t204, t157);
t145 = sin(t147);
t146 = cos(t147);
t131 = -t145 * t204 + t146 * t157;
t128 = 0.1e1 / t131;
t165 = qJ(5) + qJ(6);
t160 = cos(t165);
t205 = t166 * t160;
t159 = sin(t165);
t207 = t159 * t167;
t142 = t157 * t205 + t207;
t138 = 0.1e1 / t142;
t129 = 0.1e1 / t131 ^ 2;
t139 = 0.1e1 / t142 ^ 2;
t225 = t148 - 0.1e1;
t163 = t166 ^ 2;
t208 = t156 * t163;
t127 = t129 * t208 + 0.1e1;
t201 = qJD(1) * t167;
t183 = t156 * t166 * t201;
t200 = qJD(3) * t157;
t212 = t146 * t158;
t117 = (-t122 * t167 + qJD(3)) * t212 + (t157 * t226 + t191) * t145;
t223 = t117 * t128 * t129;
t224 = (-t208 * t223 + (-t158 * t163 * t200 + t183) * t129) / t127 ^ 2;
t161 = qJD(5) + qJD(6);
t186 = qJD(1) * t157 + t161;
t199 = qJD(3) * t166;
t175 = t158 * t199 + t186 * t167;
t187 = t157 * t161 + qJD(1);
t180 = t160 * t187;
t123 = t175 * t159 + t166 * t180;
t203 = t167 * t160;
t206 = t166 * t159;
t141 = t157 * t206 - t203;
t137 = t141 ^ 2;
t136 = t137 * t139 + 0.1e1;
t215 = t139 * t141;
t181 = t159 * t187;
t124 = t175 * t160 - t166 * t181;
t220 = t124 * t138 * t139;
t222 = (t123 * t215 - t137 * t220) / t136 ^ 2;
t221 = t122 * t158;
t219 = t129 * t158;
t218 = t129 * t166;
t210 = t152 * t158;
t178 = qJD(3) * (-t152 * t227 - t210);
t217 = (-t153 * t183 + t164 * t178) / t150 ^ 2;
t216 = t138 * t159;
t214 = t141 * t160;
t213 = t145 * t167;
t211 = t152 * t156;
t197 = -0.2e1 * t223;
t196 = 0.2e1 * t222;
t195 = t158 * t224;
t194 = t158 * t218;
t193 = t158 * t217;
t192 = t148 * t211;
t190 = 0.1e1 + t209;
t189 = 0.2e1 * t141 * t220;
t188 = t217 * t228;
t185 = t167 * t192;
t184 = t225 * t158 * t145;
t182 = t190 * t166;
t179 = t139 * t214 - t216;
t177 = t179 * t166;
t176 = t158 * t198 - t186 * t166;
t144 = t157 * t203 - t206;
t143 = t157 * t207 + t205;
t134 = 0.1e1 / t136;
t133 = t190 * t167 * t148;
t125 = 0.1e1 / t127;
t121 = (-t146 * t185 - t184) * t166;
t120 = t157 * t213 + t212 + (-t145 * t157 - t146 * t204) * t133;
t118 = -t190 * t188 + (-qJD(1) * t182 + t178 * t228) * t148;
t115 = -0.2e1 * t222 + 0.2e1 * (t123 * t134 * t139 + (-t134 * t220 - t139 * t222) * t141) * t141;
t1 = [-0.2e1 * t166 * t152 * t193 + (-qJD(3) * t182 + t201 * t210) * t148, 0, t118, 0, 0, 0; (0.2e1 * t128 * t195 + (t128 * t200 + (qJD(1) * t121 + t117) * t219) * t125) * t167 + (-0.2e1 * t129 * t195 * t121 + (((t122 * t185 + t225 * t200 + 0.2e1 * t193) * t145 + (t188 * t211 + t221 + (-t221 + (0.2e1 * t158 + t227) * t198) * t148) * t146) * t194 + (-t129 * t200 + t158 * t197) * t121 + (t128 + ((t163 - t164) * t146 * t192 - t167 * t184) * t129) * t158 * qJD(1)) * t125) * t166, 0, 0.2e1 * (-t120 * t219 - t128 * t157) * t166 * t224 + ((t128 * t201 + (-qJD(3) * t120 - t117) * t218) * t157 + (t128 * t199 + (-t118 * t146 * t167 + t226 * t145 + (-qJD(3) * t145 + t122 * t213 + t146 * t202) * t133) * t194 + (t129 * t201 + t166 * t197) * t120 + ((-t118 - t202) * t145 + ((t133 * t167 - 0.1e1) * qJD(3) + (-t133 + t167) * t122) * t146) * t157 * t218) * t158) * t125, 0, 0, 0; (-t138 * t143 + t144 * t215) * t196 + (t144 * t189 + t167 * t138 * t180 + t176 * t216 + (t167 * t141 * t181 - t144 * t123 - t143 * t124 - t176 * t214) * t139) * t134, 0, t158 * t177 * t196 + (t177 * t200 + (-t179 * t201 + ((t138 * t161 + t189) * t160 + (-t123 * t160 + (t141 * t161 - t124) * t159) * t139) * t166) * t158) * t134, 0, t115, t115;];
JaD_rot  = t1;
