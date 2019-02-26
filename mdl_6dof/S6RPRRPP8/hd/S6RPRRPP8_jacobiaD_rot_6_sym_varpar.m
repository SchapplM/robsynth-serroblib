% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_jacobiaD_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:00:21
% EndTime: 2019-02-26 21:00:22
% DurationCPUTime: 1.00s
% Computational Cost: add. (1786->113), mult. (5988->254), div. (1156->14), fcn. (7606->9), ass. (0->108)
t148 = sin(qJ(3));
t150 = cos(qJ(4));
t152 = cos(qJ(1));
t199 = t150 * t152;
t147 = sin(qJ(4));
t149 = sin(qJ(1));
t201 = t149 * t147;
t132 = t148 * t199 - t201;
t151 = cos(qJ(3));
t198 = t151 * t150;
t123 = atan2(t132, t198);
t117 = sin(t123);
t118 = cos(t123);
t128 = t132 ^ 2;
t142 = 0.1e1 / t150 ^ 2;
t145 = 0.1e1 / t151 ^ 2;
t206 = t142 * t145;
t125 = t128 * t206 + 0.1e1;
t121 = 0.1e1 / t125;
t141 = 0.1e1 / t150;
t204 = t145 * t148;
t180 = t141 * t204;
t166 = t132 * t180 + t152;
t108 = t166 * t121;
t197 = t108 - t152;
t226 = t197 * t117 * t151 + t118 * t148;
t144 = 0.1e1 / t151;
t191 = qJD(4) * t148;
t167 = (-qJD(1) - t191) * t150;
t169 = qJD(1) * t148 + qJD(4);
t193 = qJD(3) * t151;
t221 = t149 * t193 + t152 * t169;
t115 = -t221 * t147 + t149 * t167;
t129 = -t148 * t201 + t199;
t126 = t129 ^ 2;
t139 = 0.1e1 / t149 ^ 2;
t209 = t139 * t145;
t124 = t126 * t209 + 0.1e1;
t138 = 0.1e1 / t149;
t140 = t138 * t139;
t146 = t144 * t145;
t194 = qJD(3) * t148;
t177 = t146 * t194;
t195 = qJD(1) * t152;
t179 = t145 * t195;
t218 = (t129 * t115 * t209 + (t139 * t177 - t140 * t179) * t126) / t124 ^ 2;
t225 = -0.2e1 * t144 * t218;
t200 = t149 * t150;
t203 = t147 * t152;
t131 = t148 * t203 + t200;
t205 = t142 * t147;
t164 = -t131 * t141 + t132 * t205;
t224 = t144 * t164;
t222 = qJD(3) * (0.2e1 * t148 ^ 2 * t146 + t144);
t214 = t117 * t132;
t112 = t118 * t198 + t214;
t109 = 0.1e1 / t112;
t110 = 0.1e1 / t112 ^ 2;
t130 = t148 * t200 + t203;
t127 = t130 ^ 2;
t107 = t110 * t127 + 0.1e1;
t174 = t147 * t191;
t196 = qJD(1) * t149;
t116 = -t147 * t196 - t149 * t174 + t221 * t150;
t216 = t110 * t130;
t190 = qJD(4) * t151;
t163 = -t147 * t190 - t150 * t194;
t182 = t132 * t206;
t175 = t152 * t193;
t114 = qJD(1) * t130 + t131 * qJD(4) - t150 * t175;
t207 = t141 * t144;
t184 = t114 * t207;
t101 = (-t163 * t182 - t184) * t121;
t162 = -t101 * t132 - t163;
t97 = (-t101 * t198 - t114) * t117 - t162 * t118;
t219 = t109 * t110 * t97;
t220 = 0.1e1 / t107 ^ 2 * (t116 * t216 - t127 * t219);
t143 = t141 * t142;
t192 = qJD(4) * t147;
t217 = (-t114 * t182 + (t143 * t145 * t192 + t142 * t177) * t128) / t125 ^ 2;
t215 = t117 * t130;
t212 = t118 * t130;
t211 = t118 * t132;
t208 = t139 * t152;
t202 = t149 * t109;
t189 = 0.2e1 * t220;
t188 = 0.2e1 * t219;
t187 = -0.2e1 * t217;
t185 = t110 * t215;
t183 = t132 * t207;
t181 = t138 * t204;
t178 = t145 * t194;
t173 = -0.2e1 * t109 * t220;
t172 = t110 * t189;
t171 = t130 * t188;
t168 = 0.2e1 * t207 * t217;
t165 = t129 * t208 + t131 * t138;
t161 = t116 * t207 - (-t142 * t144 * t192 - t141 * t178) * t130;
t119 = 0.1e1 / t124;
t113 = t152 * t167 + (t149 * t169 - t175) * t147;
t105 = 0.1e1 / t107;
t104 = t121 * t224;
t100 = (-t117 + (-t118 * t183 + t117) * t121) * t130;
t99 = t108 * t211 - t226 * t150;
t98 = -t118 * t151 * t147 - t117 * t131 + (-t117 * t198 + t211) * t104;
t96 = t166 * t187 + (-t114 * t180 - t196 + (t141 * t222 + t174 * t206) * t132) * t121;
t94 = t187 * t224 + (t164 * t178 + (-t114 * t205 + t113 * t141 + (-t131 * t205 + (0.2e1 * t143 * t147 ^ 2 + t141) * t132) * qJD(4)) * t144) * t121;
t1 = [-t121 * t161 + t130 * t168, 0, t96, t94, 0, 0; t132 * t173 + (-t114 * t109 + (-t100 * t116 - t132 * t97) * t110) * t105 + (t100 * t172 + (t100 * t188 + (-t116 * t121 + t116 - (t101 * t121 * t183 + t187) * t130) * t110 * t117 + (-(t132 * t168 - t101) * t216 + (-(t101 + t184) * t130 + t161 * t132) * t110 * t121) * t118) * t105) * t130, 0, t99 * t130 * t172 + (t99 * t171 + (-(t211 * t96 + (-t101 * t214 - t114 * t118) * t108) * t130 - t99 * t116) * t110 + (-t151 * t202 - t226 * t216) * t192) * t105 + (t149 * t173 * t151 + ((-qJD(3) * t202 - (qJD(3) * t197 + t101) * t185) * t148 + (t109 * t195 + (-t149 * t97 - (-t96 - t196) * t215 - (-t101 * t197 - qJD(3)) * t212) * t110) * t151) * t105) * t150 (-t109 * t129 + t216 * t98) * t189 + (t98 * t171 + t115 * t109 - (t113 + (t101 * t147 - t150 * t94) * t151 + t162 * t104) * t185 + (-t98 * t116 - t129 * t97 - (t147 * t194 - t150 * t190 - t104 * t114 + t132 * t94 + (-t104 * t198 - t131) * t101) * t212) * t110) * t105, 0, 0; t165 * t225 + (t165 * t178 + (t115 * t208 - t113 * t138 + (-t131 * t208 + (-0.2e1 * t140 * t152 ^ 2 - t138) * t129) * qJD(1)) * t144) * t119, 0, 0.2e1 * (t129 * t181 - t147) * t218 + (-t115 * t181 + qJD(4) * t150 + (t139 * t148 * t179 - t138 * t222) * t129) * t119, t130 * t138 * t225 + (t116 * t138 * t144 + (-t139 * t144 * t195 + t138 * t178) * t130) * t119, 0, 0;];
JaD_rot  = t1;
