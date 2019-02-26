% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function JaD_rot = S6RPRRPP8_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_jacobiaD_rot_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:00:20
% EndTime: 2019-02-26 21:00:21
% DurationCPUTime: 0.98s
% Computational Cost: add. (1786->114), mult. (5988->255), div. (1156->14), fcn. (7606->9), ass. (0->107)
t150 = sin(qJ(3));
t151 = sin(qJ(1));
t152 = cos(qJ(4));
t202 = t151 * t152;
t149 = sin(qJ(4));
t154 = cos(qJ(1));
t205 = t149 * t154;
t131 = t150 * t205 + t202;
t153 = cos(qJ(3));
t201 = t153 * t149;
t123 = atan2(t131, t201);
t117 = sin(t123);
t118 = cos(t123);
t128 = t131 ^ 2;
t140 = 0.1e1 / t149 ^ 2;
t147 = 0.1e1 / t153 ^ 2;
t210 = t140 * t147;
t125 = t128 * t210 + 0.1e1;
t121 = 0.1e1 / t125;
t139 = 0.1e1 / t149;
t206 = t147 * t150;
t184 = t139 * t206;
t168 = t131 * t184 + t154;
t108 = t168 * t121;
t199 = t108 - t154;
t226 = -t199 * t153 * t117 - t118 * t150;
t146 = 0.1e1 / t153;
t171 = qJD(1) * t150 + qJD(4);
t169 = t171 * t154;
t172 = qJD(4) * t150 + qJD(1);
t195 = qJD(3) * t153;
t116 = t152 * t169 + (-t172 * t149 + t152 * t195) * t151;
t130 = t150 * t202 + t205;
t127 = t130 ^ 2;
t144 = 0.1e1 / t151 ^ 2;
t208 = t144 * t147;
t124 = t127 * t208 + 0.1e1;
t143 = 0.1e1 / t151;
t145 = t143 * t144;
t148 = t146 * t147;
t196 = qJD(3) * t150;
t180 = t148 * t196;
t197 = qJD(1) * t154;
t182 = t147 * t197;
t220 = (t130 * t116 * t208 + (t144 * t180 - t145 * t182) * t127) / t124 ^ 2;
t225 = -0.2e1 * t146 * t220;
t223 = qJD(3) * (0.2e1 * t150 ^ 2 * t148 + t146);
t216 = t117 * t131;
t112 = t118 * t201 + t216;
t109 = 0.1e1 / t112;
t110 = 0.1e1 / t112 ^ 2;
t200 = t154 * t152;
t203 = t151 * t149;
t129 = t150 * t203 - t200;
t126 = t129 ^ 2;
t107 = t110 * t126 + 0.1e1;
t115 = t172 * t202 + (t151 * t195 + t169) * t149;
t218 = t110 * t129;
t193 = qJD(4) * t153;
t165 = -t149 * t196 + t152 * t193;
t185 = t131 * t210;
t194 = qJD(4) * t152;
t178 = t150 * t194;
t179 = t154 * t195;
t113 = -t149 * t179 - t152 * t197 - t154 * t178 + t171 * t203;
t211 = t139 * t146;
t187 = t113 * t211;
t101 = (-t165 * t185 - t187) * t121;
t164 = -t101 * t131 - t165;
t97 = (-t101 * t201 - t113) * t117 - t164 * t118;
t221 = t109 * t110 * t97;
t222 = 0.1e1 / t107 ^ 2 * (t115 * t218 - t126 * t221);
t141 = t139 * t140;
t219 = (-t113 * t185 + (-t141 * t147 * t194 + t140 * t180) * t128) / t125 ^ 2;
t217 = t117 * t129;
t214 = t118 * t129;
t213 = t118 * t131;
t209 = t140 * t152;
t207 = t144 * t154;
t204 = t151 * t109;
t198 = qJD(1) * t151;
t192 = 0.2e1 * t222;
t191 = 0.2e1 * t221;
t190 = -0.2e1 * t219;
t188 = t110 * t217;
t186 = t131 * t211;
t183 = t143 * t206;
t181 = t147 * t196;
t177 = -0.2e1 * t109 * t222;
t176 = t110 * t192;
t175 = t129 * t191;
t174 = 0.2e1 * t146 * t219;
t170 = t139 * t174;
t132 = t150 * t200 - t203;
t167 = t130 * t207 - t132 * t143;
t166 = t131 * t209 - t132 * t139;
t163 = t115 * t211 - (t140 * t146 * t194 - t139 * t181) * t129;
t119 = 0.1e1 / t124;
t114 = qJD(1) * t130 + qJD(4) * t131 - t152 * t179;
t105 = 0.1e1 / t107;
t104 = t166 * t146 * t121;
t100 = (-t117 + (-t118 * t186 + t117) * t121) * t129;
t99 = t108 * t213 + t226 * t149;
t98 = t118 * t153 * t152 + t117 * t132 - (-t117 * t201 + t213) * t104;
t96 = t168 * t190 + (-t113 * t184 - t198 + (t139 * t223 - t178 * t210) * t131) * t121;
t94 = t166 * t174 + (-t166 * t181 + (t113 * t209 - t114 * t139 + (-t132 * t209 + (0.2e1 * t141 * t152 ^ 2 + t139) * t131) * qJD(4)) * t146) * t121;
t1 = [-t163 * t121 + t129 * t170, 0, t96, t94, 0, 0; t131 * t177 + (-t113 * t109 + (-t100 * t115 - t131 * t97) * t110) * t105 + (t100 * t176 + (t100 * t191 + (-t115 * t121 + t115 - (t101 * t121 * t186 + t190) * t129) * t110 * t117 + (-(t131 * t170 - t101) * t218 + (-(t101 + t187) * t129 + t163 * t131) * t110 * t121) * t118) * t105) * t129, 0, t99 * t129 * t176 + (t99 * t175 + (-(t96 * t213 + (-t101 * t216 - t113 * t118) * t108) * t129 - t99 * t115) * t110 + (t153 * t204 - t226 * t218) * t194) * t105 + (t151 * t177 * t153 + ((-qJD(3) * t204 - (t199 * qJD(3) + t101) * t188) * t150 + (t109 * t197 + (-t151 * t97 - (-t96 - t198) * t217 - (-t199 * t101 - qJD(3)) * t214) * t110) * t153) * t105) * t149 (-t109 * t130 + t98 * t218) * t192 + (t98 * t175 + t116 * t109 - (-t114 + (-t101 * t152 - t149 * t94) * t153 - t164 * t104) * t188 + (-t98 * t115 - t130 * t97 - (-t152 * t196 - t149 * t193 + t104 * t113 + t131 * t94 + (t104 * t201 + t132) * t101) * t214) * t110) * t105, 0, 0; t167 * t225 + (t167 * t181 + (t116 * t207 + t114 * t143 + (t132 * t207 + (-0.2e1 * t145 * t154 ^ 2 - t143) * t130) * qJD(1)) * t146) * t119, 0, 0.2e1 * (t130 * t183 + t152) * t220 + (-t116 * t183 + qJD(4) * t149 + (t144 * t150 * t182 - t143 * t223) * t130) * t119, t129 * t143 * t225 + (t115 * t143 * t146 + (-t144 * t146 * t197 + t143 * t181) * t129) * t119, 0, 0;];
JaD_rot  = t1;
