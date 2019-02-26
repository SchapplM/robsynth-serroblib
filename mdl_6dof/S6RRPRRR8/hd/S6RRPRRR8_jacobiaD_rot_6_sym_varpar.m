% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR8
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
% Datum: 2019-02-26 21:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:58:16
% EndTime: 2019-02-26 21:58:17
% DurationCPUTime: 0.72s
% Computational Cost: add. (2878->96), mult. (2949->204), div. (516->12), fcn. (3430->9), ass. (0->95)
t162 = sin(qJ(2));
t156 = t162 ^ 2;
t164 = cos(qJ(2));
t159 = 0.1e1 / t164 ^ 2;
t209 = t156 * t159;
t163 = sin(qJ(1));
t227 = 0.2e1 * t163;
t226 = t162 * t209;
t153 = pkin(11) + qJ(4) + qJ(5) + qJ(6);
t152 = cos(t153);
t165 = cos(qJ(1));
t201 = t164 * t165;
t151 = sin(t153);
t205 = t163 * t151;
t141 = t152 * t201 + t205;
t203 = t163 * t162;
t146 = atan2(-t203, -t164);
t145 = cos(t146);
t144 = sin(t146);
t190 = t144 * t203;
t131 = -t145 * t164 - t190;
t128 = 0.1e1 / t131;
t135 = 0.1e1 / t141;
t158 = 0.1e1 / t164;
t129 = 0.1e1 / t131 ^ 2;
t136 = 0.1e1 / t141 ^ 2;
t225 = -0.2e1 * t162;
t157 = t163 ^ 2;
t150 = t157 * t209 + 0.1e1;
t148 = 0.1e1 / t150;
t224 = t148 - 0.1e1;
t199 = qJD(1) * t165;
t187 = t162 * t199;
t197 = qJD(2) * t164;
t198 = qJD(2) * t163;
t124 = (-(-t163 * t197 - t187) * t158 + t198 * t209) * t148;
t212 = t145 * t162;
t115 = (-t124 * t163 + qJD(2)) * t212 + (-t187 + (t124 - t198) * t164) * t144;
t223 = t115 * t128 * t129;
t154 = qJD(4) + qJD(5) + qJD(6);
t180 = -qJD(1) * t164 + t154;
t181 = t154 * t164 - qJD(1);
t196 = qJD(2) * t165;
t186 = t162 * t196;
t211 = t151 * t165;
t120 = -t181 * t211 + (t180 * t163 - t186) * t152;
t222 = t120 * t135 * t136;
t202 = t163 * t164;
t175 = t151 * t202 + t152 * t165;
t119 = t175 * qJD(1) - t141 * t154 + t151 * t186;
t204 = t163 * t152;
t140 = t151 * t201 - t204;
t134 = t140 ^ 2;
t123 = t134 * t136 + 0.1e1;
t214 = t136 * t140;
t221 = 0.1e1 / t123 ^ 2 * (-t119 * t214 - t134 * t222);
t220 = t124 * t144;
t219 = t124 * t162;
t218 = t129 * t162;
t207 = t158 * t162;
t174 = qJD(2) * (t158 * t226 + t207);
t178 = t156 * t163 * t199;
t217 = (t157 * t174 + t159 * t178) / t150 ^ 2;
t185 = 0.1e1 + t209;
t133 = t185 * t163 * t148;
t216 = t133 * t163;
t215 = t135 * t151;
t213 = t140 * t152;
t210 = t156 * t158;
t161 = t165 ^ 2;
t208 = t156 * t161;
t206 = t162 * t165;
t200 = qJD(1) * t163;
t127 = t129 * t208 + 0.1e1;
t195 = 0.2e1 * (-t208 * t223 + (t161 * t162 * t197 - t178) * t129) / t127 ^ 2;
t194 = 0.2e1 * t223;
t193 = -0.2e1 * t221;
t192 = t140 * t222;
t191 = t129 * t206;
t189 = t148 * t210;
t184 = t162 * t195;
t183 = t217 * t225;
t182 = t217 * t227;
t179 = t163 * t189;
t177 = t185 * t165;
t176 = t136 * t213 - t215;
t173 = t162 * t198 + t180 * t165;
t139 = -t152 * t202 + t211;
t125 = 0.1e1 / t127;
t121 = 0.1e1 / t123;
t118 = (t224 * t162 * t144 - t145 * t179) * t165;
t117 = -t144 * t202 + t212 + (t144 * t164 - t145 * t203) * t133;
t116 = -t185 * t182 + (qJD(1) * t177 + t174 * t227) * t148;
t112 = t193 + 0.2e1 * (-t119 * t136 * t121 + (-t121 * t222 - t136 * t221) * t140) * t140;
t1 = [t165 * t158 * t183 + (qJD(2) * t177 - t200 * t207) * t148, t116, 0, 0, 0, 0; (t128 * t184 + (-t128 * t197 + (qJD(1) * t118 + t115) * t218) * t125) * t163 + (t129 * t184 * t118 + (-((t124 * t179 + t224 * t197 + t183) * t144 + (t182 * t210 - t219 + (t219 + (t225 - t226) * t198) * t148) * t145) * t191 + (-t129 * t197 + t162 * t194) * t118 + (-t128 + ((-t157 + t161) * t145 * t189 + t224 * t190) * t129) * t162 * qJD(1)) * t125) * t165 (t117 * t218 - t128 * t164) * t165 * t195 + ((-t128 * t200 + (-qJD(2) * t117 - t115) * t165 * t129) * t164 + (-t128 * t196 - (-t116 * t145 * t163 + t144 * t198 + t216 * t220 - t220 + (-qJD(2) * t144 - t145 * t199) * t133) * t191 + (t129 * t200 + t165 * t194) * t117 - ((t116 - t199) * t144 + ((0.1e1 - t216) * qJD(2) + (t133 - t163) * t124) * t145) * t129 * t201) * t162) * t125, 0, 0, 0, 0; 0.2e1 * (t135 * t175 + t139 * t214) * t221 + (0.2e1 * t139 * t192 - t181 * t135 * t204 + t173 * t215 + (-t181 * t140 * t205 + t139 * t119 + t120 * t175 - t173 * t213) * t136) * t121, t176 * t193 * t206 + (t176 * t164 * t196 + (-t176 * t200 + ((-t135 * t154 - 0.2e1 * t192) * t152 + (-t119 * t152 + (-t140 * t154 + t120) * t151) * t136) * t165) * t162) * t121, 0, t112, t112, t112;];
JaD_rot  = t1;
