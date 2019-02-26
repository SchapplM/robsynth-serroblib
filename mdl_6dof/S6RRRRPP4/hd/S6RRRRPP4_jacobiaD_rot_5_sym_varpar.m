% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:26:59
% EndTime: 2019-02-26 22:27:00
% DurationCPUTime: 0.76s
% Computational Cost: add. (2009->96), mult. (2734->204), div. (498->12), fcn. (3199->9), ass. (0->95)
t157 = sin(qJ(2));
t151 = t157 ^ 2;
t159 = cos(qJ(2));
t154 = 0.1e1 / t159 ^ 2;
t204 = t151 * t154;
t158 = sin(qJ(1));
t222 = 0.2e1 * t158;
t221 = t157 * t204;
t148 = qJ(3) + qJ(4) + pkin(10);
t147 = cos(t148);
t160 = cos(qJ(1));
t196 = t159 * t160;
t146 = sin(t148);
t200 = t158 * t146;
t136 = t147 * t196 + t200;
t198 = t158 * t157;
t141 = atan2(-t198, -t159);
t139 = cos(t141);
t138 = sin(t141);
t185 = t138 * t198;
t126 = -t139 * t159 - t185;
t123 = 0.1e1 / t126;
t130 = 0.1e1 / t136;
t153 = 0.1e1 / t159;
t124 = 0.1e1 / t126 ^ 2;
t131 = 0.1e1 / t136 ^ 2;
t220 = -0.2e1 * t157;
t152 = t158 ^ 2;
t144 = t152 * t204 + 0.1e1;
t142 = 0.1e1 / t144;
t219 = t142 - 0.1e1;
t149 = qJD(3) + qJD(4);
t197 = t158 * t159;
t170 = t146 * t197 + t147 * t160;
t191 = qJD(2) * t160;
t181 = t157 * t191;
t114 = t170 * qJD(1) - t136 * t149 + t146 * t181;
t199 = t158 * t147;
t135 = t146 * t196 - t199;
t129 = t135 ^ 2;
t119 = t129 * t131 + 0.1e1;
t209 = t131 * t135;
t175 = -qJD(1) * t159 + t149;
t176 = t149 * t159 - qJD(1);
t206 = t146 * t160;
t115 = -t176 * t206 + (t175 * t158 - t181) * t147;
t216 = t115 * t130 * t131;
t218 = (-t114 * t209 - t129 * t216) / t119 ^ 2;
t194 = qJD(1) * t160;
t182 = t157 * t194;
t192 = qJD(2) * t159;
t193 = qJD(2) * t158;
t116 = (-(-t158 * t192 - t182) * t153 + t193 * t204) * t142;
t207 = t139 * t157;
t110 = (-t116 * t158 + qJD(2)) * t207 + (-t182 + (t116 - t193) * t159) * t138;
t217 = t110 * t123 * t124;
t215 = t116 * t138;
t214 = t116 * t157;
t213 = t124 * t157;
t202 = t153 * t157;
t169 = qJD(2) * (t153 * t221 + t202);
t173 = t151 * t158 * t194;
t212 = (t152 * t169 + t154 * t173) / t144 ^ 2;
t180 = 0.1e1 + t204;
t128 = t180 * t158 * t142;
t211 = t128 * t158;
t210 = t130 * t146;
t208 = t135 * t147;
t205 = t151 * t153;
t156 = t160 ^ 2;
t203 = t151 * t156;
t201 = t157 * t160;
t195 = qJD(1) * t158;
t122 = t124 * t203 + 0.1e1;
t190 = 0.2e1 * (-t203 * t217 + (t156 * t157 * t192 - t173) * t124) / t122 ^ 2;
t189 = -0.2e1 * t218;
t188 = 0.2e1 * t217;
t187 = t135 * t216;
t186 = t124 * t201;
t184 = t142 * t205;
t179 = t157 * t190;
t178 = t212 * t220;
t177 = t212 * t222;
t174 = t158 * t184;
t172 = t180 * t160;
t171 = t131 * t208 - t210;
t168 = t157 * t193 + t175 * t160;
t134 = -t147 * t197 + t206;
t120 = 0.1e1 / t122;
t117 = 0.1e1 / t119;
t113 = (t219 * t157 * t138 - t139 * t174) * t160;
t112 = -t138 * t197 + t207 + (t138 * t159 - t139 * t198) * t128;
t111 = -t180 * t177 + (qJD(1) * t172 + t169 * t222) * t142;
t107 = t189 + 0.2e1 * (-t114 * t117 * t131 + (-t117 * t216 - t131 * t218) * t135) * t135;
t1 = [t153 * t160 * t178 + (qJD(2) * t172 - t195 * t202) * t142, t111, 0, 0, 0, 0; (t123 * t179 + (-t123 * t192 + (qJD(1) * t113 + t110) * t213) * t120) * t158 + (t124 * t179 * t113 + (-((t116 * t174 + t219 * t192 + t178) * t138 + (t177 * t205 - t214 + (t214 + (t220 - t221) * t193) * t142) * t139) * t186 + (-t124 * t192 + t157 * t188) * t113 + (-t123 + ((-t152 + t156) * t139 * t184 + t219 * t185) * t124) * t157 * qJD(1)) * t120) * t160 (t112 * t213 - t123 * t159) * t160 * t190 + ((-t123 * t195 + (-qJD(2) * t112 - t110) * t160 * t124) * t159 + (-t123 * t191 - (-t111 * t139 * t158 + t138 * t193 + t211 * t215 - t215 + (-qJD(2) * t138 - t139 * t194) * t128) * t186 + (t124 * t195 + t160 * t188) * t112 - ((t111 - t194) * t138 + ((0.1e1 - t211) * qJD(2) + (t128 - t158) * t116) * t139) * t124 * t196) * t157) * t120, 0, 0, 0, 0; 0.2e1 * (t130 * t170 + t134 * t209) * t218 + (0.2e1 * t134 * t187 - t176 * t130 * t199 + t168 * t210 + (-t176 * t135 * t200 + t134 * t114 + t115 * t170 - t168 * t208) * t131) * t117, t171 * t189 * t201 + (t171 * t159 * t191 + (-t171 * t195 + ((-t130 * t149 - 0.2e1 * t187) * t147 + (-t114 * t147 + (-t135 * t149 + t115) * t146) * t131) * t160) * t157) * t117, t107, t107, 0, 0;];
JaD_rot  = t1;
