% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:40
% EndTime: 2019-02-26 21:37:41
% DurationCPUTime: 0.66s
% Computational Cost: add. (5265->87), mult. (3414->191), div. (723->12), fcn. (4027->9), ass. (0->91)
t154 = sin(qJ(1));
t209 = 0.2e1 * t154;
t148 = qJ(2) + pkin(10) + qJ(4);
t147 = cos(t148);
t153 = cos(pkin(11));
t183 = t154 * t153;
t152 = sin(pkin(11));
t155 = cos(qJ(1));
t187 = t152 * t155;
t135 = t147 * t187 - t183;
t129 = t135 ^ 2;
t184 = t154 * t152;
t186 = t153 * t155;
t136 = t147 * t186 + t184;
t131 = 0.1e1 / t136 ^ 2;
t124 = t129 * t131 + 0.1e1;
t133 = -t147 * t184 - t186;
t146 = sin(t148);
t149 = qJD(2) + qJD(4);
t188 = t149 * t155;
t173 = t146 * t188;
t125 = qJD(1) * t133 - t152 * t173;
t197 = t131 * t135;
t134 = -t147 * t183 + t187;
t126 = qJD(1) * t134 - t153 * t173;
t130 = 0.1e1 / t136;
t201 = t126 * t130 * t131;
t208 = (t125 * t197 - t129 * t201) / t124 ^ 2;
t150 = t154 ^ 2;
t142 = t146 ^ 2;
t144 = 0.1e1 / t147 ^ 2;
t195 = t142 * t144;
t140 = t150 * t195 + 0.1e1;
t138 = 0.1e1 / t140;
t143 = 0.1e1 / t147;
t181 = qJD(1) * t155;
t172 = t146 * t181;
t189 = t149 * t154;
t175 = t144 * t189;
t112 = (-(-t147 * t189 - t172) * t143 + t142 * t175) * t138;
t207 = t112 - t189;
t185 = t154 * t146;
t137 = atan2(-t185, -t147);
t128 = cos(t137);
t127 = sin(t137);
t176 = t127 * t185;
t120 = -t128 * t147 - t176;
t117 = 0.1e1 / t120;
t118 = 0.1e1 / t120 ^ 2;
t206 = t138 - 0.1e1;
t199 = t128 * t146;
t107 = (-t112 * t154 + t149) * t199 + (t207 * t147 - t172) * t127;
t205 = t107 * t117 * t118;
t141 = t146 * t142;
t192 = t143 * t146;
t163 = t149 * (t141 * t143 * t144 + t192);
t193 = t142 * t154;
t166 = t181 * t193;
t204 = (t144 * t166 + t150 * t163) / t140 ^ 2;
t203 = t118 * t146;
t202 = t118 * t155;
t200 = t127 * t154;
t198 = t130 * t152;
t196 = t142 * t143;
t151 = t155 ^ 2;
t194 = t142 * t151;
t191 = t146 * t155;
t190 = t147 * t149;
t182 = qJD(1) * t154;
t115 = t118 * t194 + 0.1e1;
t180 = 0.2e1 * (-t194 * t205 + (t146 * t151 * t190 - t166) * t118) / t115 ^ 2;
t179 = 0.2e1 * t205;
t178 = t118 * t191;
t177 = t135 * t201;
t174 = t149 * t185;
t171 = 0.1e1 + t195;
t170 = t146 * t180;
t169 = -0.2e1 * t146 * t204;
t168 = t204 * t209;
t167 = t128 * t138 * t196;
t165 = t171 * t155;
t164 = t153 * t197 - t198;
t122 = 0.1e1 / t124;
t121 = t171 * t154 * t138;
t113 = 0.1e1 / t115;
t111 = (t206 * t146 * t127 - t154 * t167) * t155;
t109 = -t147 * t200 + t199 + (t127 * t147 - t128 * t185) * t121;
t108 = -t171 * t168 + (qJD(1) * t165 + t163 * t209) * t138;
t105 = -0.2e1 * t164 * t191 * t208 + (t164 * t147 * t188 + (-0.2e1 * t177 * t186 + t182 * t198 + (t126 * t187 + (t125 * t155 - t135 * t182) * t153) * t131) * t146) * t122;
t104 = (t109 * t203 - t117 * t147) * t155 * t180 + ((-t117 * t182 + (-t109 * t149 - t107) * t202) * t147 + (-t117 * t188 - (-t108 * t128 * t154 - t207 * t127 + (t112 * t200 - t127 * t149 - t128 * t181) * t121) * t178 + (t118 * t182 + t155 * t179) * t109 - ((t108 - t181) * t127 + ((-t121 * t154 + 0.1e1) * t149 + (t121 - t154) * t112) * t128) * t147 * t202) * t146) * t113;
t1 = [t143 * t155 * t169 + (t149 * t165 - t182 * t192) * t138, t108, 0, t108, 0, 0; (t117 * t170 + (-t117 * t190 + (qJD(1) * t111 + t107) * t203) * t113) * t154 + (t118 * t170 * t111 + (-((t169 - t190 + (t112 * t143 * t193 + t190) * t138) * t127 + (t168 * t196 - t112 * t146 + (-t141 * t175 + (t112 - 0.2e1 * t189) * t146) * t138) * t128) * t178 + (-t118 * t190 + t146 * t179) * t111 + (-t117 + ((-t150 + t151) * t167 + t206 * t176) * t118) * t146 * qJD(1)) * t113) * t155, t104, 0, t104, 0, 0; 0.2e1 * (-t130 * t133 + t134 * t197) * t208 + ((-qJD(1) * t135 + t152 * t174) * t130 + 0.2e1 * t134 * t177 + (-t133 * t126 - (-qJD(1) * t136 + t153 * t174) * t135 - t134 * t125) * t131) * t122, t105, 0, t105, 0, 0;];
JaD_rot  = t1;
