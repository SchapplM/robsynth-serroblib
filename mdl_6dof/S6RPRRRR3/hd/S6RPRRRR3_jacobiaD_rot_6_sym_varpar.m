% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JaD_rot = S6RPRRRR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:16:00
% EndTime: 2019-02-26 21:16:01
% DurationCPUTime: 0.80s
% Computational Cost: add. (3504->97), mult. (2949->204), div. (516->12), fcn. (3430->9), ass. (0->96)
t162 = sin(qJ(3));
t158 = t162 ^ 2;
t163 = cos(qJ(3));
t160 = 0.1e1 / t163 ^ 2;
t201 = t158 * t160;
t156 = qJ(1) + pkin(11);
t152 = sin(t156);
t225 = 0.2e1 * t152;
t224 = t162 * t201;
t155 = qJ(4) + qJ(5) + qJ(6);
t149 = cos(t155);
t153 = cos(t156);
t203 = t153 * t163;
t148 = sin(t155);
t208 = t152 * t148;
t138 = t149 * t203 + t208;
t206 = t152 * t162;
t143 = atan2(-t206, -t163);
t142 = cos(t143);
t141 = sin(t143);
t189 = t141 * t206;
t128 = -t142 * t163 - t189;
t125 = 0.1e1 / t128;
t132 = 0.1e1 / t138;
t159 = 0.1e1 / t163;
t126 = 0.1e1 / t128 ^ 2;
t133 = 0.1e1 / t138 ^ 2;
t223 = -0.2e1 * t162;
t150 = t152 ^ 2;
t146 = t150 * t201 + 0.1e1;
t144 = 0.1e1 / t146;
t222 = t144 - 0.1e1;
t154 = qJD(4) + qJD(5) + qJD(6);
t205 = t152 * t163;
t174 = t148 * t205 + t153 * t149;
t196 = qJD(3) * t162;
t185 = t153 * t196;
t116 = t174 * qJD(1) - t138 * t154 + t148 * t185;
t207 = t152 * t149;
t137 = t148 * t203 - t207;
t131 = t137 ^ 2;
t121 = t131 * t133 + 0.1e1;
t212 = t133 * t137;
t179 = -qJD(1) * t163 + t154;
t180 = t154 * t163 - qJD(1);
t204 = t153 * t148;
t117 = -t180 * t204 + (t179 * t152 - t185) * t149;
t219 = t117 * t132 * t133;
t221 = (-t116 * t212 - t131 * t219) / t121 ^ 2;
t198 = qJD(1) * t162;
t186 = t153 * t198;
t195 = qJD(3) * t163;
t197 = qJD(3) * t152;
t118 = (-(-t152 * t195 - t186) * t159 + t197 * t201) * t144;
t210 = t142 * t162;
t112 = (-t118 * t152 + qJD(3)) * t210 + (-t186 + (t118 - t197) * t163) * t141;
t220 = t112 * t125 * t126;
t218 = t118 * t141;
t217 = t118 * t162;
t216 = t126 * t162;
t173 = qJD(3) * (t162 + t224) * t159;
t199 = qJD(1) * t153;
t177 = t152 * t158 * t199;
t215 = (t150 * t173 + t160 * t177) / t146 ^ 2;
t184 = 0.1e1 + t201;
t130 = t184 * t152 * t144;
t214 = t130 * t152;
t213 = t132 * t148;
t211 = t137 * t149;
t151 = t153 ^ 2;
t209 = t151 * t158;
t202 = t158 * t159;
t200 = qJD(1) * t152;
t124 = t126 * t209 + 0.1e1;
t194 = 0.2e1 * (-t209 * t220 + (t151 * t162 * t195 - t177) * t126) / t124 ^ 2;
t193 = 0.2e1 * t221;
t192 = 0.2e1 * t220;
t191 = t137 * t219;
t190 = t153 * t216;
t188 = t144 * t202;
t183 = t162 * t194;
t182 = t215 * t225;
t181 = t215 * t223;
t178 = t152 * t188;
t176 = t184 * t153;
t175 = t133 * t211 - t213;
t172 = t175 * t162;
t171 = t152 * t196 + t179 * t153;
t136 = -t149 * t205 + t204;
t122 = 0.1e1 / t124;
t119 = 0.1e1 / t121;
t115 = (t222 * t162 * t141 - t142 * t178) * t153;
t114 = -t141 * t205 + t210 + (t141 * t163 - t142 * t206) * t130;
t113 = -t184 * t182 + (qJD(1) * t176 + t173 * t225) * t144;
t109 = -0.2e1 * t221 + 0.2e1 * (-t116 * t119 * t133 + (-t119 * t219 - t133 * t221) * t137) * t137;
t1 = [t153 * t159 * t181 + (-t152 * t159 * t198 + qJD(3) * t176) * t144, 0, t113, 0, 0, 0; (t125 * t183 + (-t125 * t195 + (qJD(1) * t115 + t112) * t216) * t122) * t152 + (t126 * t183 * t115 + (-((t118 * t178 + t222 * t195 + t181) * t141 + (t182 * t202 - t217 + (t217 + (t223 - t224) * t197) * t144) * t142) * t190 + (-t126 * t195 + t162 * t192) * t115 + (-t125 + ((-t150 + t151) * t142 * t188 + t222 * t189) * t126) * t198) * t122) * t153, 0 (t114 * t216 - t125 * t163) * t153 * t194 + ((-t125 * t200 + (-qJD(3) * t114 - t112) * t153 * t126) * t163 + (-t153 * qJD(3) * t125 - (-t113 * t142 * t152 + t141 * t197 + t214 * t218 - t218 + (-qJD(3) * t141 - t142 * t199) * t130) * t190 + (t126 * t200 + t153 * t192) * t114 - ((t113 - t199) * t141 + ((0.1e1 - t214) * qJD(3) + (t130 - t152) * t118) * t142) * t126 * t203) * t162) * t122, 0, 0, 0; (t132 * t174 + t136 * t212) * t193 + (0.2e1 * t136 * t191 - t180 * t132 * t207 + t171 * t213 + (-t180 * t137 * t208 + t136 * t116 + t117 * t174 - t171 * t211) * t133) * t119, 0, -t153 * t172 * t193 + (-t172 * t200 + (t175 * t195 + ((-t132 * t154 - 0.2e1 * t191) * t149 + (-t116 * t149 + (-t137 * t154 + t117) * t148) * t133) * t162) * t153) * t119, t109, t109, t109;];
JaD_rot  = t1;
