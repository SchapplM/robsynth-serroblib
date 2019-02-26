% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:01:26
% EndTime: 2019-02-26 21:01:26
% DurationCPUTime: 0.74s
% Computational Cost: add. (3074->97), mult. (2734->204), div. (498->12), fcn. (3199->9), ass. (0->96)
t158 = sin(qJ(3));
t154 = t158 ^ 2;
t159 = cos(qJ(3));
t156 = 0.1e1 / t159 ^ 2;
t197 = t154 * t156;
t152 = qJ(1) + pkin(10);
t148 = sin(t152);
t221 = 0.2e1 * t148;
t220 = t158 * t197;
t150 = qJ(4) + pkin(11) + qJ(6);
t145 = cos(t150);
t149 = cos(t152);
t199 = t149 * t159;
t144 = sin(t150);
t204 = t148 * t144;
t134 = t145 * t199 + t204;
t202 = t148 * t158;
t139 = atan2(-t202, -t159);
t137 = cos(t139);
t136 = sin(t139);
t185 = t136 * t202;
t124 = -t137 * t159 - t185;
t121 = 0.1e1 / t124;
t128 = 0.1e1 / t134;
t155 = 0.1e1 / t159;
t122 = 0.1e1 / t124 ^ 2;
t129 = 0.1e1 / t134 ^ 2;
t219 = -0.2e1 * t158;
t146 = t148 ^ 2;
t142 = t146 * t197 + 0.1e1;
t140 = 0.1e1 / t142;
t218 = t140 - 0.1e1;
t151 = qJD(4) + qJD(6);
t201 = t148 * t159;
t170 = t144 * t201 + t149 * t145;
t192 = qJD(3) * t158;
t181 = t149 * t192;
t112 = t170 * qJD(1) - t134 * t151 + t144 * t181;
t203 = t148 * t145;
t133 = t144 * t199 - t203;
t127 = t133 ^ 2;
t117 = t127 * t129 + 0.1e1;
t208 = t129 * t133;
t175 = -qJD(1) * t159 + t151;
t176 = t151 * t159 - qJD(1);
t200 = t149 * t144;
t113 = -t176 * t200 + (t175 * t148 - t181) * t145;
t215 = t113 * t128 * t129;
t217 = (-t112 * t208 - t127 * t215) / t117 ^ 2;
t194 = qJD(1) * t158;
t182 = t149 * t194;
t191 = qJD(3) * t159;
t193 = qJD(3) * t148;
t114 = (-(-t148 * t191 - t182) * t155 + t193 * t197) * t140;
t206 = t137 * t158;
t108 = (-t114 * t148 + qJD(3)) * t206 + (-t182 + (t114 - t193) * t159) * t136;
t216 = t108 * t121 * t122;
t214 = t114 * t136;
t213 = t114 * t158;
t212 = t122 * t158;
t169 = qJD(3) * (t158 + t220) * t155;
t195 = qJD(1) * t149;
t173 = t148 * t154 * t195;
t211 = (t146 * t169 + t156 * t173) / t142 ^ 2;
t180 = 0.1e1 + t197;
t126 = t180 * t148 * t140;
t210 = t126 * t148;
t209 = t128 * t144;
t207 = t133 * t145;
t147 = t149 ^ 2;
t205 = t147 * t154;
t198 = t154 * t155;
t196 = qJD(1) * t148;
t120 = t122 * t205 + 0.1e1;
t190 = 0.2e1 * (-t205 * t216 + (t147 * t158 * t191 - t173) * t122) / t120 ^ 2;
t189 = 0.2e1 * t217;
t188 = 0.2e1 * t216;
t187 = t133 * t215;
t186 = t149 * t212;
t184 = t140 * t198;
t179 = t158 * t190;
t178 = t211 * t221;
t177 = t211 * t219;
t174 = t148 * t184;
t172 = t180 * t149;
t171 = t129 * t207 - t209;
t168 = t171 * t158;
t167 = t148 * t192 + t175 * t149;
t132 = -t145 * t201 + t200;
t118 = 0.1e1 / t120;
t115 = 0.1e1 / t117;
t111 = (t218 * t158 * t136 - t137 * t174) * t149;
t110 = -t136 * t201 + t206 + (t136 * t159 - t137 * t202) * t126;
t109 = -t180 * t178 + (qJD(1) * t172 + t169 * t221) * t140;
t105 = -0.2e1 * t217 + 0.2e1 * (-t112 * t115 * t129 + (-t115 * t215 - t129 * t217) * t133) * t133;
t1 = [t149 * t155 * t177 + (-t148 * t155 * t194 + qJD(3) * t172) * t140, 0, t109, 0, 0, 0; (t121 * t179 + (-t121 * t191 + (qJD(1) * t111 + t108) * t212) * t118) * t148 + (t122 * t179 * t111 + (-((t114 * t174 + t218 * t191 + t177) * t136 + (t178 * t198 - t213 + (t213 + (t219 - t220) * t193) * t140) * t137) * t186 + (-t122 * t191 + t158 * t188) * t111 + (-t121 + ((-t146 + t147) * t137 * t184 + t218 * t185) * t122) * t194) * t118) * t149, 0 (t110 * t212 - t121 * t159) * t149 * t190 + ((-t121 * t196 + (-qJD(3) * t110 - t108) * t149 * t122) * t159 + (-t149 * qJD(3) * t121 - (-t109 * t137 * t148 + t136 * t193 + t210 * t214 - t214 + (-qJD(3) * t136 - t137 * t195) * t126) * t186 + (t122 * t196 + t149 * t188) * t110 - ((t109 - t195) * t136 + ((0.1e1 - t210) * qJD(3) + (t126 - t148) * t114) * t137) * t122 * t199) * t158) * t118, 0, 0, 0; (t128 * t170 + t132 * t208) * t189 + (0.2e1 * t132 * t187 - t176 * t128 * t203 + t167 * t209 + (-t176 * t133 * t204 + t132 * t112 + t113 * t170 - t167 * t207) * t129) * t115, 0, -t149 * t168 * t189 + (-t168 * t196 + (t171 * t191 + ((-t128 * t151 - 0.2e1 * t187) * t145 + (-t112 * t145 + (-t133 * t151 + t113) * t144) * t129) * t158) * t149) * t115, t105, 0, t105;];
JaD_rot  = t1;
