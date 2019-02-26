% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:57:25
% EndTime: 2019-02-26 20:57:27
% DurationCPUTime: 1.15s
% Computational Cost: add. (3939->115), mult. (5988->257), div. (1156->14), fcn. (7606->9), ass. (0->107)
t141 = qJ(1) + pkin(9);
t139 = sin(t141);
t140 = cos(t141);
t149 = sin(qJ(4));
t151 = cos(qJ(4));
t152 = cos(qJ(3));
t203 = t151 * t152;
t126 = t139 * t203 - t140 * t149;
t150 = sin(qJ(3));
t204 = t150 * t151;
t120 = atan2(-t126, t204);
t116 = sin(t120);
t117 = cos(t120);
t122 = t126 ^ 2;
t143 = 0.1e1 / t150 ^ 2;
t146 = 0.1e1 / t151 ^ 2;
t208 = t143 * t146;
t121 = t122 * t208 + 0.1e1;
t118 = 0.1e1 / t121;
t145 = 0.1e1 / t151;
t207 = t143 * t152;
t182 = t145 * t207;
t166 = t126 * t182 + t139;
t104 = t166 * t118;
t226 = t104 - t139;
t227 = t226 * t116 * t150 - t117 * t152;
t194 = qJD(4) * t152;
t176 = t149 * t194;
t198 = qJD(3) * t150;
t225 = t151 * t198 + t176;
t142 = 0.1e1 / t150;
t144 = t142 * t143;
t224 = qJD(3) * (0.2e1 * t144 * t152 ^ 2 + t142);
t199 = qJD(1) * t152;
t181 = t139 * t199;
t195 = qJD(4) * t151;
t200 = qJD(1) * t140;
t110 = -t139 * t195 + t225 * t140 - t149 * t200 + t151 * t181;
t129 = t139 * t149 + t140 * t203;
t197 = qJD(3) * t152;
t180 = t143 * t197;
t196 = qJD(4) * t149;
t209 = t142 * t145;
t223 = (-t142 * t146 * t196 + t145 * t180) * t129 + t110 * t209;
t222 = (-qJD(1) + t194) * t151 - t149 * t198;
t218 = t116 * t126;
t108 = t117 * t204 - t218;
t105 = 0.1e1 / t108;
t136 = 0.1e1 / t140;
t106 = 0.1e1 / t108 ^ 2;
t137 = 0.1e1 / t140 ^ 2;
t163 = -t150 * t196 + t151 * t197;
t184 = t126 * t208;
t112 = t129 * qJD(1) - t225 * t139 - t140 * t195;
t186 = t112 * t209;
t97 = (t163 * t184 - t186) * t118;
t161 = -t126 * t97 + t163;
t168 = -t97 * t204 - t112;
t93 = t168 * t116 + t161 * t117;
t221 = t105 * t106 * t93;
t147 = t145 * t146;
t179 = t144 * t197;
t220 = 0.1e1 / t121 ^ 2 * (t112 * t184 + (t143 * t147 * t196 - t146 * t179) * t122);
t219 = t106 * t129;
t217 = t116 * t129;
t215 = t117 * t126;
t214 = t117 * t129;
t212 = t137 * t139;
t211 = t137 * t143;
t210 = t140 * t105;
t206 = t146 * t149;
t205 = t149 * t152;
t201 = qJD(1) * t139;
t124 = t129 ^ 2;
t103 = t124 * t106 + 0.1e1;
t193 = 0.2e1 / t103 ^ 2 * (-t110 * t219 - t124 * t221);
t192 = 0.2e1 * t221;
t167 = (-qJD(4) + t199) * t149;
t109 = t139 * t167 - t222 * t140;
t128 = t139 * t151 - t140 * t205;
t123 = t128 ^ 2;
t115 = t123 * t211 + 0.1e1;
t138 = t136 * t137;
t191 = 0.2e1 / t115 ^ 2 * (t128 * t109 * t211 + (t138 * t143 * t201 - t137 * t179) * t123);
t190 = -0.2e1 * t220;
t189 = t142 * t220;
t188 = t106 * t217;
t185 = t126 * t209;
t183 = t136 * t207;
t175 = t105 * t193;
t174 = t106 * t193;
t173 = t129 * t192;
t172 = t142 * t191;
t170 = t145 * t189;
t125 = t139 * t205 + t140 * t151;
t165 = -t125 * t145 + t126 * t206;
t164 = -t125 * t136 - t128 * t212;
t113 = 0.1e1 / t115;
t111 = t222 * t139 + t140 * t167;
t101 = 0.1e1 / t103;
t100 = t165 * t142 * t118;
t96 = (-t116 + (t117 * t185 + t116) * t118) * t129;
t95 = -t104 * t215 - t227 * t151;
t94 = -t117 * t150 * t149 + t116 * t125 - (-t116 * t204 - t215) * t100;
t92 = t166 * t190 + (t112 * t182 + t200 + (-t145 * t224 + t176 * t208) * t126) * t118;
t90 = 0.2e1 * t165 * t189 + (t165 * t180 + (-t112 * t206 + t111 * t145 + (t125 * t206 + (-0.2e1 * t147 * t149 ^ 2 - t145) * t126) * qJD(4)) * t142) * t118;
t1 = [t223 * t118 + 0.2e1 * t129 * t170, 0, t92, t90, 0, 0; t126 * t175 + (-t112 * t105 + (t110 * t96 + t126 * t93) * t106) * t101 + (t96 * t174 + (t96 * t192 + (t110 * t118 - t110 - (-t118 * t97 * t185 + t190) * t129) * t106 * t116 + (-(-0.2e1 * t126 * t170 - t97) * t219 + (-(t97 + t186) * t129 + t223 * t126) * t106 * t118) * t117) * t101) * t129, 0, t95 * t129 * t174 + (t95 * t173 + (-(-t92 * t215 + (-t112 * t117 + t218 * t97) * t104) * t129 + t95 * t110) * t106 + (t150 * t210 - t227 * t219) * t196) * t101 + (t140 * t175 * t150 + ((-qJD(3) * t210 - (-qJD(3) * t226 - t97) * t188) * t152 + (t105 * t201 + (t140 * t93 - (-t92 + t200) * t217 - (-t226 * t97 - qJD(3)) * t214) * t106) * t150) * t101) * t151 (-t105 * t128 + t94 * t219) * t193 + (t94 * t173 + t109 * t105 - (t111 + (t149 * t97 - t151 * t90) * t150 + t161 * t100) * t188 + (t94 * t110 - t128 * t93 - (-t168 * t100 + t125 * t97 - t126 * t90 - t149 * t197 - t150 * t195) * t214) * t106) * t101, 0, 0; t164 * t172 + (t164 * t180 + (t109 * t212 + t111 * t136 + (t125 * t212 + (0.2e1 * t138 * t139 ^ 2 + t136) * t128) * qJD(1)) * t142) * t113, 0 (t128 * t183 - t149) * t191 + (-t109 * t183 + t195 + (t136 * t224 - t181 * t211) * t128) * t113, t129 * t136 * t172 + (t110 * t136 * t142 + (-t137 * t142 * t201 + t136 * t180) * t129) * t113, 0, 0;];
JaD_rot  = t1;
