% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:44:51
% EndTime: 2019-02-26 20:44:52
% DurationCPUTime: 0.98s
% Computational Cost: add. (7031->125), mult. (6168->273), div. (1114->15), fcn. (7752->9), ass. (0->115)
t161 = pkin(10) + qJ(5);
t159 = cos(t161);
t158 = sin(t161);
t214 = qJ(1) + pkin(9);
t196 = sin(t214);
t189 = t196 * t158;
t160 = cos(t214);
t167 = cos(qJ(3));
t222 = t160 * t167;
t143 = t159 * t222 + t189;
t137 = 0.1e1 / t143 ^ 2;
t157 = t160 ^ 2;
t166 = sin(qJ(3));
t162 = t166 ^ 2;
t225 = t157 * t162;
t205 = t137 * t225;
t132 = 0.1e1 + t205;
t186 = qJD(1) * t196;
t219 = qJD(3) * t160;
t201 = t166 * t219;
t175 = t167 * t186 + t201;
t185 = t196 * qJD(5);
t224 = t160 * t158;
t122 = (-qJD(5) * t167 + qJD(1)) * t224 + (t185 - t175) * t159;
t136 = 0.1e1 / t143;
t235 = t122 * t136 * t137;
t191 = t225 * t235;
t218 = qJD(3) * t167;
t199 = t166 * t218;
t242 = (-t191 + (-t160 * t162 * t186 + t157 * t199) * t137) / t132 ^ 2;
t223 = t160 * t166;
t139 = t160 * t159 + t167 * t189;
t182 = t158 * t185;
t216 = qJD(5) * t160;
t197 = t159 * t216;
t121 = t139 * qJD(1) + t158 * t201 - t167 * t197 - t182;
t188 = t196 * t159;
t142 = t158 * t222 - t188;
t154 = 0.1e1 / t158;
t155 = 0.1e1 / t158 ^ 2;
t163 = 0.1e1 / t166;
t164 = 0.1e1 / t166 ^ 2;
t200 = t164 * t218;
t217 = qJD(5) * t159;
t228 = t154 * t163;
t241 = (t155 * t163 * t217 + t154 * t200) * t142 + t121 * t228;
t221 = t166 * t158;
t131 = atan2(-t139, t221);
t126 = cos(t131);
t125 = sin(t131);
t234 = t125 * t139;
t120 = t126 * t221 - t234;
t117 = 0.1e1 / t120;
t118 = 0.1e1 / t120 ^ 2;
t240 = 0.2e1 * t142;
t134 = t139 ^ 2;
t226 = t155 * t164;
t133 = t134 * t226 + 0.1e1;
t129 = 0.1e1 / t133;
t215 = qJD(5) * t166;
t180 = t158 * t218 + t159 * t215;
t203 = t139 * t226;
t190 = t166 * t196;
t183 = qJD(3) * t190;
t184 = t159 * t186;
t220 = qJD(1) * t160;
t123 = t159 * t185 * t167 - t184 + (t220 * t167 - t183 - t216) * t158;
t206 = t123 * t228;
t109 = (t180 * t203 - t206) * t129;
t176 = -t109 * t139 + t180;
t105 = (-t109 * t221 - t123) * t125 + t176 * t126;
t119 = t117 * t118;
t239 = t105 * t119;
t156 = t154 * t155;
t165 = t163 / t162;
t198 = t164 * t217;
t238 = (t123 * t203 + (-t155 * t165 * t218 - t156 * t198) * t134) / t133 ^ 2;
t237 = t118 * t142;
t236 = t121 * t118;
t233 = t125 * t142;
t232 = t125 * t166;
t231 = t126 * t139;
t230 = t126 * t142;
t229 = t126 * t167;
t227 = t155 * t159;
t135 = t142 ^ 2;
t115 = t118 * t135 + 0.1e1;
t213 = 0.2e1 * (-t135 * t239 - t142 * t236) / t115 ^ 2;
t212 = 0.2e1 * t242;
t211 = -0.2e1 * t238;
t210 = t119 * t240;
t209 = t163 * t238;
t208 = t118 * t233;
t204 = t139 * t228;
t202 = t154 * t164 * t167;
t195 = t117 * t213;
t194 = t118 * t213;
t193 = t223 * t240;
t192 = t154 * t209;
t179 = t139 * t202 + t196;
t116 = t179 * t129;
t187 = t196 - t116;
t141 = t167 * t188 - t224;
t181 = t139 * t227 - t141 * t154;
t178 = t137 * t141 * t160 - t196 * t136;
t127 = 0.1e1 / t132;
t124 = t143 * qJD(1) - t159 * t183 - t167 * t182 - t197;
t113 = 0.1e1 / t115;
t112 = t181 * t163 * t129;
t108 = (-t125 + (t126 * t204 + t125) * t129) * t142;
t107 = -t116 * t231 + (t187 * t232 + t229) * t158;
t106 = t126 * t159 * t166 - t125 * t141 + (-t125 * t221 - t231) * t112;
t104 = t179 * t211 + (t123 * t202 + t220 + (-t155 * t167 * t198 + (-0.2e1 * t165 * t167 ^ 2 - t163) * t154 * qJD(3)) * t139) * t129;
t102 = -0.2e1 * t181 * t209 + (-t181 * t200 + (t123 * t227 - t124 * t154 + (t141 * t227 + (-0.2e1 * t156 * t159 ^ 2 - t154) * t139) * qJD(5)) * t163) * t129;
t1 = [t241 * t129 + t192 * t240, 0, t104, 0, t102, 0; t139 * t195 + (-t123 * t117 + (t105 * t139 + t108 * t121) * t118) * t113 + (t108 * t194 + (0.2e1 * t108 * t239 + (t121 * t129 - t121 - (-t109 * t129 * t204 + t211) * t142) * t118 * t125 + (-(-0.2e1 * t139 * t192 - t109) * t237 + (-(t109 + t206) * t142 + t241 * t139) * t118 * t129) * t126) * t113) * t142, 0, t107 * t142 * t194 + (-(-t104 * t231 + (t109 * t234 - t123 * t126) * t116) * t237 + (t105 * t210 + t236) * t107 + (-t117 * t223 - (-t116 * t232 + t125 * t190 + t229) * t237) * t217) * t113 + (t195 * t223 + ((-t117 * t219 - (t187 * qJD(3) - t109) * t208) * t167 + (t117 * t186 + (t160 * t105 - (-t104 + t220) * t233 - (t187 * t109 - qJD(3)) * t230) * t118) * t166) * t113) * t158, 0 (t106 * t237 - t117 * t143) * t213 + (t106 * t236 + t122 * t117 + (t106 * t210 - t143 * t118) * t105 - (t159 * t218 - t158 * t215 - t102 * t139 - t112 * t123 + (-t112 * t221 - t141) * t109) * t118 * t230 - (-t124 + (-t102 * t158 - t109 * t159) * t166 - t176 * t112) * t208) * t113, 0; t178 * t166 * t212 + (-t178 * t218 + ((qJD(1) * t136 + 0.2e1 * t141 * t235) * t160 + (-t196 * t122 - t124 * t160 + t141 * t186) * t137) * t166) * t127, 0 (t136 * t222 + t159 * t205) * t212 + (0.2e1 * t159 * t191 + t175 * t136 + ((t122 * t167 + 0.2e1 * t162 * t184) * t160 + (qJD(5) * t158 * t162 - 0.2e1 * t159 * t199) * t157) * t137) * t127, 0, t137 * t193 * t242 + (t193 * t235 + (t121 * t223 + (-t160 * t218 + t166 * t186) * t142) * t137) * t127, 0;];
JaD_rot  = t1;
