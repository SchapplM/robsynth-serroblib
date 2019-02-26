% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRP6_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP6_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:42:35
% EndTime: 2019-02-26 22:42:35
% DurationCPUTime: 0.71s
% Computational Cost: add. (2348->96), mult. (2949->204), div. (516->12), fcn. (3430->9), ass. (0->95)
t159 = sin(qJ(2));
t153 = t159 ^ 2;
t161 = cos(qJ(2));
t156 = 0.1e1 / t161 ^ 2;
t206 = t153 * t156;
t160 = sin(qJ(1));
t224 = 0.2e1 * t160;
t223 = t159 * t206;
t151 = qJ(3) + qJ(4) + qJ(5);
t149 = cos(t151);
t162 = cos(qJ(1));
t198 = t161 * t162;
t148 = sin(t151);
t202 = t160 * t148;
t138 = t149 * t198 + t202;
t200 = t160 * t159;
t143 = atan2(-t200, -t161);
t142 = cos(t143);
t141 = sin(t143);
t187 = t141 * t200;
t128 = -t142 * t161 - t187;
t125 = 0.1e1 / t128;
t132 = 0.1e1 / t138;
t155 = 0.1e1 / t161;
t126 = 0.1e1 / t128 ^ 2;
t133 = 0.1e1 / t138 ^ 2;
t222 = -0.2e1 * t159;
t154 = t160 ^ 2;
t146 = t154 * t206 + 0.1e1;
t144 = 0.1e1 / t146;
t221 = t144 - 0.1e1;
t150 = qJD(3) + qJD(4) + qJD(5);
t199 = t160 * t161;
t172 = t148 * t199 + t149 * t162;
t193 = qJD(2) * t162;
t183 = t159 * t193;
t116 = t172 * qJD(1) - t138 * t150 + t148 * t183;
t201 = t160 * t149;
t137 = t148 * t198 - t201;
t131 = t137 ^ 2;
t121 = t131 * t133 + 0.1e1;
t211 = t133 * t137;
t177 = -qJD(1) * t161 + t150;
t178 = t150 * t161 - qJD(1);
t208 = t148 * t162;
t117 = -t178 * t208 + (t177 * t160 - t183) * t149;
t218 = t117 * t132 * t133;
t220 = (-t116 * t211 - t131 * t218) / t121 ^ 2;
t196 = qJD(1) * t162;
t184 = t159 * t196;
t194 = qJD(2) * t161;
t195 = qJD(2) * t160;
t118 = (-(-t160 * t194 - t184) * t155 + t195 * t206) * t144;
t209 = t142 * t159;
t112 = (-t118 * t160 + qJD(2)) * t209 + (-t184 + (t118 - t195) * t161) * t141;
t219 = t112 * t125 * t126;
t217 = t118 * t141;
t216 = t118 * t159;
t215 = t126 * t159;
t204 = t155 * t159;
t171 = qJD(2) * (t155 * t223 + t204);
t175 = t153 * t160 * t196;
t214 = (t154 * t171 + t156 * t175) / t146 ^ 2;
t182 = 0.1e1 + t206;
t130 = t182 * t160 * t144;
t213 = t130 * t160;
t212 = t132 * t148;
t210 = t137 * t149;
t207 = t153 * t155;
t158 = t162 ^ 2;
t205 = t153 * t158;
t203 = t159 * t162;
t197 = qJD(1) * t160;
t124 = t126 * t205 + 0.1e1;
t192 = 0.2e1 * (-t205 * t219 + (t158 * t159 * t194 - t175) * t126) / t124 ^ 2;
t191 = -0.2e1 * t220;
t190 = 0.2e1 * t219;
t189 = t137 * t218;
t188 = t126 * t203;
t186 = t144 * t207;
t181 = t159 * t192;
t180 = t214 * t222;
t179 = t214 * t224;
t176 = t160 * t186;
t174 = t182 * t162;
t173 = t133 * t210 - t212;
t170 = t159 * t195 + t177 * t162;
t136 = -t149 * t199 + t208;
t122 = 0.1e1 / t124;
t119 = 0.1e1 / t121;
t115 = (t221 * t159 * t141 - t142 * t176) * t162;
t114 = -t141 * t199 + t209 + (t141 * t161 - t142 * t200) * t130;
t113 = -t182 * t179 + (qJD(1) * t174 + t171 * t224) * t144;
t109 = t191 + 0.2e1 * (-t116 * t119 * t133 + (-t119 * t218 - t133 * t220) * t137) * t137;
t1 = [t155 * t162 * t180 + (qJD(2) * t174 - t197 * t204) * t144, t113, 0, 0, 0, 0; (t125 * t181 + (-t125 * t194 + (qJD(1) * t115 + t112) * t215) * t122) * t160 + (t126 * t181 * t115 + (-((t118 * t176 + t221 * t194 + t180) * t141 + (t179 * t207 - t216 + (t216 + (t222 - t223) * t195) * t144) * t142) * t188 + (-t126 * t194 + t159 * t190) * t115 + (-t125 + ((-t154 + t158) * t142 * t186 + t221 * t187) * t126) * t159 * qJD(1)) * t122) * t162 (t114 * t215 - t125 * t161) * t162 * t192 + ((-t125 * t197 + (-qJD(2) * t114 - t112) * t162 * t126) * t161 + (-t125 * t193 - (-t113 * t142 * t160 + t141 * t195 + t213 * t217 - t217 + (-qJD(2) * t141 - t142 * t196) * t130) * t188 + (t126 * t197 + t162 * t190) * t114 - ((t113 - t196) * t141 + ((0.1e1 - t213) * qJD(2) + (t130 - t160) * t118) * t142) * t126 * t198) * t159) * t122, 0, 0, 0, 0; 0.2e1 * (t132 * t172 + t136 * t211) * t220 + (0.2e1 * t136 * t189 - t178 * t132 * t201 + t170 * t212 + (-t178 * t137 * t202 + t136 * t116 + t117 * t172 - t170 * t210) * t133) * t119, t173 * t191 * t203 + (t173 * t161 * t193 + (-t173 * t197 + ((-t132 * t150 - 0.2e1 * t189) * t149 + (-t116 * t149 + (-t137 * t150 + t117) * t148) * t133) * t162) * t159) * t119, t109, t109, t109, 0;];
JaD_rot  = t1;
