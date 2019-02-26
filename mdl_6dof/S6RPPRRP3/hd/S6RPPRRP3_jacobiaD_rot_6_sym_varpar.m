% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRP3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:31:18
% EndTime: 2019-02-26 20:31:19
% DurationCPUTime: 0.97s
% Computational Cost: add. (4034->122), mult. (6168->269), div. (1114->15), fcn. (7752->9), ass. (0->112)
t142 = qJ(1) + pkin(9);
t141 = sin(t142);
t153 = cos(qJ(4));
t207 = t141 * t153;
t151 = sin(qJ(4));
t150 = sin(qJ(5));
t224 = cos(t142);
t186 = t224 * t150;
t152 = cos(qJ(5));
t208 = t141 * t152;
t132 = t151 * t186 + t208;
t203 = t153 * t150;
t123 = atan2(t132, t203);
t120 = cos(t123);
t119 = sin(t123);
t216 = t119 * t132;
t111 = t120 * t203 + t216;
t108 = 0.1e1 / t111;
t131 = t151 * t208 + t186;
t126 = 0.1e1 / t131;
t143 = 0.1e1 / t150;
t147 = 0.1e1 / t153;
t109 = 0.1e1 / t111 ^ 2;
t127 = 0.1e1 / t131 ^ 2;
t144 = 0.1e1 / t150 ^ 2;
t185 = t224 * t152;
t209 = t141 * t151;
t130 = t150 * t209 - t185;
t125 = t130 ^ 2;
t106 = t125 * t109 + 0.1e1;
t180 = t224 * qJD(1);
t169 = t151 * t180;
t163 = t224 * qJD(5) + t169;
t174 = qJD(5) * t151 + qJD(1);
t200 = qJD(4) * t153;
t114 = t163 * t150 + (t150 * t200 + t174 * t152) * t141;
t219 = t109 * t130;
t129 = t132 ^ 2;
t148 = 0.1e1 / t153 ^ 2;
t205 = t144 * t148;
t124 = t129 * t205 + 0.1e1;
t121 = 0.1e1 / t124;
t198 = qJD(5) * t153;
t201 = qJD(4) * t151;
t164 = -t150 * t201 + t152 * t198;
t190 = t132 * t205;
t187 = t153 * t224;
t170 = qJD(4) * t187;
t171 = t152 * t180;
t172 = t151 * t185;
t210 = t141 * t150;
t112 = -qJD(5) * t172 - t150 * t170 - t171 + (qJD(1) * t151 + qJD(5)) * t210;
t206 = t143 * t147;
t192 = t112 * t206;
t100 = (-t164 * t190 - t192) * t121;
t162 = -t100 * t132 - t164;
t96 = (-t100 * t203 - t112) * t119 - t162 * t120;
t222 = t108 * t109 * t96;
t223 = 0.1e1 / t106 ^ 2 * (t114 * t219 - t125 * t222);
t140 = t141 ^ 2;
t146 = t153 ^ 2;
t211 = t140 * t146;
t189 = t127 * t211;
t118 = 0.1e1 + t189;
t183 = t152 * t200;
t115 = t163 * t152 + (-t174 * t150 + t183) * t141;
t218 = t115 * t126 * t127;
t173 = t211 * t218;
t221 = (-t173 + (-t140 * t151 * t200 + t141 * t146 * t180) * t127) / t118 ^ 2;
t145 = t143 * t144;
t149 = t147 / t146;
t199 = qJD(5) * t152;
t182 = t148 * t199;
t220 = (-t112 * t190 + (t144 * t149 * t201 - t145 * t182) * t129) / t124 ^ 2;
t217 = t119 * t130;
t215 = t119 * t153;
t214 = t120 * t130;
t213 = t120 * t132;
t212 = t120 * t151;
t204 = t144 * t152;
t202 = qJD(1) * t141;
t197 = 0.2e1 * t223;
t196 = 0.2e1 * t222;
t195 = 0.2e1 * t221;
t194 = -0.2e1 * t220;
t193 = t109 * t217;
t191 = t132 * t206;
t188 = t143 * t148 * t151;
t184 = t148 * t201;
t166 = t132 * t188 + t224;
t107 = t166 * t121;
t181 = t224 - t107;
t179 = -0.2e1 * t108 * t223;
t178 = t109 * t197;
t177 = t130 * t196;
t176 = 0.2e1 * t147 * t220;
t175 = -0.2e1 * t130 * t207;
t168 = t143 * t176;
t133 = t172 - t210;
t167 = t132 * t204 - t133 * t143;
t165 = t127 * t133 * t141 - t224 * t126;
t161 = t114 * t206 - (t144 * t147 * t199 - t143 * t184) * t130;
t116 = 0.1e1 / t118;
t113 = t131 * qJD(1) + t132 * qJD(5) - t152 * t170;
t104 = 0.1e1 / t106;
t103 = t167 * t147 * t121;
t99 = (-t119 + (-t120 * t191 + t119) * t121) * t130;
t98 = t107 * t213 + (t181 * t215 - t212) * t150;
t97 = t120 * t153 * t152 + t119 * t133 - (-t119 * t203 + t213) * t103;
t95 = t166 * t194 + (-t112 * t188 - t202 + (-t144 * t151 * t182 + (0.2e1 * t149 * t151 ^ 2 + t147) * t143 * qJD(4)) * t132) * t121;
t93 = t167 * t176 + (-t167 * t184 + (t112 * t204 - t113 * t143 + (-t133 * t204 + (0.2e1 * t145 * t152 ^ 2 + t143) * t132) * qJD(5)) * t147) * t121;
t1 = [-t161 * t121 + t130 * t168, 0, 0, t95, t93, 0; t132 * t179 + (-t112 * t108 + (-t114 * t99 - t132 * t96) * t109) * t104 + (t99 * t178 + (t99 * t196 + (-t114 * t121 + t114 - (t100 * t121 * t191 + t194) * t130) * t109 * t119 + (-(t132 * t168 - t100) * t219 + (-(t100 + t192) * t130 + t161 * t132) * t109 * t121) * t120) * t104) * t130, 0, 0, t98 * t130 * t178 + (t98 * t177 + (-(t95 * t213 + (-t100 * t216 - t112 * t120) * t107) * t130 - t98 * t114) * t109 + (t108 * t207 - (-t107 * t215 + t119 * t187 - t212) * t219) * t199) * t104 + (t179 * t207 + ((-t141 * qJD(4) * t108 - (-t181 * qJD(4) + t100) * t193) * t151 + (t108 * t180 + (-t141 * t96 - (-t95 - t202) * t217 - (t181 * t100 - qJD(4)) * t214) * t109) * t153) * t104) * t150 (-t108 * t131 + t97 * t219) * t197 + (t97 * t177 + t115 * t108 - (-t113 + (-t100 * t152 - t150 * t93) * t153 - t162 * t103) * t193 + (-t97 * t114 - t131 * t96 - (-t152 * t201 - t150 * t198 + t103 * t112 + t132 * t93 + (t103 * t203 + t133) * t100) * t214) * t109) * t104, 0; t165 * t153 * t195 + (t165 * t201 + ((-qJD(1) * t126 + 0.2e1 * t133 * t218) * t141 + (t113 * t141 - t224 * t115 - t133 * t180) * t127) * t153) * t116, 0, 0 (t126 * t209 + t152 * t189) * t195 + (0.2e1 * t152 * t173 + (-t141 * t200 - t169) * t126 + ((t115 * t151 - 0.2e1 * t146 * t171) * t141 + (qJD(5) * t146 * t150 + 0.2e1 * t151 * t183) * t140) * t127) * t116, t127 * t175 * t221 + (t175 * t218 + (t114 * t207 + (-t141 * t201 + t153 * t180) * t130) * t127) * t116, 0;];
JaD_rot  = t1;
