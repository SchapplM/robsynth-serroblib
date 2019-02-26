% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:16:00
% EndTime: 2019-02-26 21:16:01
% DurationCPUTime: 0.80s
% Computational Cost: add. (2635->97), mult. (2734->204), div. (498->12), fcn. (3199->9), ass. (0->95)
t158 = sin(qJ(3));
t153 = t158 ^ 2;
t159 = cos(qJ(3));
t155 = 0.1e1 / t159 ^ 2;
t198 = t153 * t155;
t151 = qJ(1) + pkin(11);
t146 = sin(t151);
t222 = 0.2e1 * t146;
t221 = t158 * t198;
t147 = cos(t151);
t157 = qJ(4) + qJ(5);
t148 = sin(t157);
t149 = cos(t157);
t200 = t149 * t159;
t134 = t146 * t148 + t147 * t200;
t150 = qJD(4) + qJD(5);
t176 = t150 * t159 - qJD(1);
t193 = qJD(3) * t158;
t220 = t176 * t148 + t149 * t193;
t203 = t146 * t158;
t138 = atan2(-t203, -t159);
t137 = cos(t138);
t136 = sin(t138);
t186 = t136 * t203;
t124 = -t137 * t159 - t186;
t121 = 0.1e1 / t124;
t128 = 0.1e1 / t134;
t154 = 0.1e1 / t159;
t122 = 0.1e1 / t124 ^ 2;
t129 = 0.1e1 / t134 ^ 2;
t219 = -0.2e1 * t158;
t144 = t146 ^ 2;
t142 = t144 * t198 + 0.1e1;
t140 = 0.1e1 / t142;
t218 = t140 - 0.1e1;
t201 = t148 * t159;
t169 = t146 * t201 + t147 * t149;
t182 = t148 * t193;
t112 = t169 * qJD(1) - t134 * t150 + t147 * t182;
t133 = -t146 * t149 + t147 * t201;
t127 = t133 ^ 2;
t120 = t127 * t129 + 0.1e1;
t208 = t129 * t133;
t175 = -qJD(1) * t159 + t150;
t171 = t175 * t149;
t113 = t146 * t171 - t220 * t147;
t215 = t113 * t128 * t129;
t217 = (-t112 * t208 - t127 * t215) / t120 ^ 2;
t195 = qJD(1) * t158;
t183 = t147 * t195;
t192 = qJD(3) * t159;
t194 = qJD(3) * t146;
t114 = (-(-t146 * t192 - t183) * t154 + t194 * t198) * t140;
t206 = t137 * t158;
t108 = (-t114 * t146 + qJD(3)) * t206 + (-t183 + (t114 - t194) * t159) * t136;
t216 = t108 * t121 * t122;
t214 = t114 * t136;
t213 = t114 * t158;
t212 = t122 * t147;
t211 = t122 * t158;
t168 = qJD(3) * (t158 + t221) * t154;
t196 = qJD(1) * t147;
t173 = t146 * t153 * t196;
t210 = (t144 * t168 + t155 * t173) / t142 ^ 2;
t180 = 0.1e1 + t198;
t126 = t180 * t146 * t140;
t209 = t126 * t146;
t207 = t136 * t159;
t145 = t147 ^ 2;
t205 = t145 * t153;
t202 = t147 * t148;
t199 = t153 * t154;
t197 = qJD(1) * t146;
t117 = t122 * t205 + 0.1e1;
t191 = 0.2e1 * (-t205 * t216 + (t145 * t158 * t192 - t173) * t122) / t117 ^ 2;
t190 = 0.2e1 * t217;
t189 = 0.2e1 * t216;
t188 = t133 * t215;
t187 = t147 * t211;
t185 = t140 * t199;
t179 = t158 * t191;
t178 = t210 * t222;
t177 = t210 * t219;
t174 = t146 * t185;
t172 = t180 * t147;
t170 = -t128 * t148 + t149 * t208;
t167 = t170 * t158;
t132 = -t146 * t200 + t202;
t118 = 0.1e1 / t120;
t115 = 0.1e1 / t117;
t111 = (t218 * t158 * t136 - t137 * t174) * t147;
t110 = -t146 * t207 + t206 + (-t137 * t203 + t207) * t126;
t109 = -t180 * t178 + (qJD(1) * t172 + t168 * t222) * t140;
t105 = -0.2e1 * t217 + 0.2e1 * (-t112 * t118 * t129 + (-t118 * t215 - t129 * t217) * t133) * t133;
t1 = [t147 * t154 * t177 + (-t146 * t154 * t195 + qJD(3) * t172) * t140, 0, t109, 0, 0, 0; (t121 * t179 + (-t121 * t192 + (qJD(1) * t111 + t108) * t211) * t115) * t146 + (t122 * t179 * t111 + (-((t114 * t174 + t218 * t192 + t177) * t136 + (t178 * t199 - t213 + (t213 + (t219 - t221) * t194) * t140) * t137) * t187 + (-t122 * t192 + t158 * t189) * t111 + (-t121 + ((-t144 + t145) * t137 * t185 + t218 * t186) * t122) * t195) * t115) * t147, 0 (t110 * t211 - t121 * t159) * t147 * t191 + ((-t121 * t197 + (-qJD(3) * t110 - t108) * t212) * t159 + (-t147 * qJD(3) * t121 - (-t109 * t137 * t146 + t136 * t194 + t209 * t214 - t214 + (-qJD(3) * t136 - t137 * t196) * t126) * t187 + (t122 * t197 + t147 * t189) * t110 - ((t109 - t196) * t136 + ((0.1e1 - t209) * qJD(3) + (t126 - t146) * t114) * t137) * t159 * t212) * t158) * t115, 0, 0, 0; (t128 * t169 + t132 * t208) * t190 + (0.2e1 * t132 * t188 + (t132 * t112 + t169 * t113 + (-t220 * t146 - t147 * t171) * t133) * t129 + (t175 * t202 + (-t176 * t149 + t182) * t146) * t128) * t118, 0, -t147 * t167 * t190 + (-t167 * t197 + (t170 * t192 + ((-t128 * t150 - 0.2e1 * t188) * t149 + (-t112 * t149 + (-t133 * t150 + t113) * t148) * t129) * t158) * t147) * t118, t105, t105, 0;];
JaD_rot  = t1;
