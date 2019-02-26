% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:45:18
% EndTime: 2019-02-26 20:45:19
% DurationCPUTime: 1.11s
% Computational Cost: add. (4034->118), mult. (6168->265), div. (1114->15), fcn. (7752->9), ass. (0->113)
t143 = qJ(1) + pkin(9);
t142 = cos(t143);
t154 = cos(qJ(3));
t208 = t142 * t154;
t228 = 0.2e1 * t208;
t141 = sin(t143);
t152 = sin(qJ(3));
t153 = cos(qJ(5));
t203 = t152 * t153;
t151 = sin(qJ(5));
t209 = t142 * t151;
t132 = t141 * t203 + t209;
t202 = t154 * t153;
t123 = atan2(t132, t202);
t119 = sin(t123);
t120 = cos(t123);
t129 = t132 ^ 2;
t145 = 0.1e1 / t153 ^ 2;
t149 = 0.1e1 / t154 ^ 2;
t206 = t145 * t149;
t124 = t129 * t206 + 0.1e1;
t121 = 0.1e1 / t124;
t144 = 0.1e1 / t153;
t181 = t144 * t149 * t152;
t168 = t132 * t181 + t141;
t107 = t168 * t121;
t201 = t107 - t141;
t227 = t201 * t154 * t119 + t120 * t152;
t204 = t151 * t152;
t131 = t141 * t153 + t142 * t204;
t127 = 0.1e1 / t131 ^ 2;
t140 = t142 ^ 2;
t147 = t154 ^ 2;
t211 = t140 * t147;
t183 = t127 * t211;
t118 = 0.1e1 + t183;
t172 = qJD(5) * t152 + qJD(1);
t196 = qJD(3) * t154;
t164 = t151 * t196 + t172 * t153;
t198 = qJD(1) * t152;
t171 = qJD(5) + t198;
t210 = t141 * t151;
t115 = t164 * t142 - t171 * t210;
t126 = 0.1e1 / t131;
t218 = t115 * t126 * t127;
t170 = t211 * t218;
t178 = t152 * t196;
t200 = qJD(1) * t141;
t180 = t147 * t200;
t226 = (-t170 + (-t140 * t178 - t142 * t180) * t127) / t118 ^ 2;
t148 = 0.1e1 / t154;
t133 = -t141 * t204 + t142 * t153;
t205 = t145 * t151;
t166 = t132 * t205 + t133 * t144;
t225 = t148 * t166;
t216 = t119 * t132;
t111 = t120 * t202 + t216;
t108 = 0.1e1 / t111;
t109 = 0.1e1 / t111 ^ 2;
t223 = -0.2e1 * t151;
t182 = t142 * t203;
t130 = -t182 + t210;
t125 = t130 ^ 2;
t106 = t125 * t109 + 0.1e1;
t177 = t153 * t196;
t114 = t132 * qJD(1) + t131 * qJD(5) - t142 * t177;
t219 = t109 * t130;
t193 = qJD(5) * t154;
t197 = qJD(3) * t152;
t165 = -t151 * t193 - t153 * t197;
t184 = t132 * t206;
t194 = qJD(5) * t153;
t112 = -qJD(1) * t182 - t141 * t177 - t142 * t194 + t172 * t210;
t207 = t144 * t148;
t186 = t112 * t207;
t100 = (-t165 * t184 - t186) * t121;
t163 = -t100 * t132 - t165;
t96 = (-t100 * t202 - t112) * t119 - t163 * t120;
t221 = t108 * t109 * t96;
t222 = 0.1e1 / t106 ^ 2 * (t114 * t219 - t125 * t221);
t146 = t144 * t145;
t150 = t148 / t147;
t195 = qJD(5) * t151;
t176 = t149 * t195;
t220 = (-t112 * t184 + (t145 * t150 * t197 + t146 * t176) * t129) / t124 ^ 2;
t217 = t119 * t130;
t214 = t120 * t130;
t213 = t120 * t132;
t199 = qJD(1) * t142;
t192 = 0.2e1 * t222;
t191 = 0.2e1 * t221;
t190 = 0.2e1 * t226;
t189 = -0.2e1 * t220;
t188 = t108 * t222;
t187 = t109 * t217;
t185 = t132 * t207;
t179 = t149 * t197;
t175 = t109 * t192;
t174 = t130 * t191;
t173 = t130 * t228;
t169 = 0.2e1 * t207 * t220;
t167 = -t127 * t133 * t142 - t126 * t141;
t162 = t114 * t207 - (-t145 * t148 * t195 - t144 * t179) * t130;
t116 = 0.1e1 / t118;
t113 = t164 * t141 + t171 * t209;
t104 = 0.1e1 / t106;
t103 = t121 * t225;
t99 = (-t119 + (-t120 * t185 + t119) * t121) * t130;
t98 = t107 * t213 - t153 * t227;
t97 = -t120 * t154 * t151 + t119 * t133 + (-t119 * t202 + t213) * t103;
t95 = t168 * t189 + (-t112 * t181 + t199 + (t145 * t152 * t176 + (0.2e1 * t150 * t152 ^ 2 + t148) * t144 * qJD(3)) * t132) * t121;
t93 = t189 * t225 + (t166 * t179 + (-t112 * t205 - t113 * t144 + (t133 * t205 + (0.2e1 * t146 * t151 ^ 2 + t144) * t132) * qJD(5)) * t148) * t121;
t1 = [-t162 * t121 + t130 * t169, 0, t95, 0, t93, 0; -0.2e1 * t132 * t188 + (-t112 * t108 + (-t114 * t99 - t132 * t96) * t109) * t104 + (t99 * t175 + (t99 * t191 + (-t114 * t121 + t114 - (t100 * t121 * t185 + t189) * t130) * t109 * t119 + (-(t132 * t169 - t100) * t219 + (-(t100 + t186) * t130 + t162 * t132) * t109 * t121) * t120) * t104) * t130, 0, t98 * t130 * t175 + (t98 * t174 + (-(t95 * t213 + (-t100 * t216 - t112 * t120) * t107) * t130 - t98 * t114) * t109 + (t108 * t208 - t219 * t227) * t195) * t104 + (t188 * t228 + ((t142 * qJD(3) * t108 - (t201 * qJD(3) + t100) * t187) * t152 + (t108 * t200 + (t142 * t96 - (-t95 + t199) * t217 - (-t201 * t100 - qJD(3)) * t214) * t109) * t154) * t104) * t153, 0 (-t108 * t131 + t97 * t219) * t192 + (t97 * t174 + t115 * t108 - (-t113 + (t100 * t151 - t153 * t93) * t154 + t163 * t103) * t187 + (-t97 * t114 - t131 * t96 - (t151 * t197 - t153 * t193 - t103 * t112 + t132 * t93 + (-t103 * t202 + t133) * t100) * t214) * t109) * t104, 0; t167 * t154 * t190 + (t167 * t197 + ((qJD(1) * t126 - 0.2e1 * t133 * t218) * t142 + (-t113 * t142 + (-qJD(1) * t133 - t115) * t141) * t127) * t154) * t116, 0 (-t126 * t142 * t152 - t151 * t183) * t190 + (t170 * t223 + (-t141 * t198 + t142 * t196) * t126 + ((-t115 * t152 + t180 * t223) * t142 + (t147 * t194 + t178 * t223) * t140) * t127) * t116, 0, t127 * t173 * t226 + (t173 * t218 + (-t114 * t208 + (t142 * t197 + t154 * t200) * t130) * t127) * t116, 0;];
JaD_rot  = t1;
