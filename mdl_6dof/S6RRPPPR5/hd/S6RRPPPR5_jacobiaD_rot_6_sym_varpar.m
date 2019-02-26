% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:18
% EndTime: 2019-02-26 21:24:19
% DurationCPUTime: 0.84s
% Computational Cost: add. (1304->107), mult. (3345->233), div. (468->12), fcn. (4023->11), ass. (0->103)
t156 = sin(qJ(2));
t147 = t156 ^ 2;
t159 = cos(qJ(2));
t150 = 0.1e1 / t159 ^ 2;
t204 = t147 * t150;
t157 = sin(qJ(1));
t222 = 0.2e1 * t157;
t221 = t156 * t204;
t155 = sin(qJ(6));
t158 = cos(qJ(6));
t160 = cos(qJ(1));
t192 = qJD(6) * t160;
t197 = qJD(1) * t157;
t220 = t155 * t197 - t158 * t192;
t219 = t155 * t192 + t158 * t197;
t153 = sin(pkin(9));
t154 = cos(pkin(9));
t198 = t159 * t160;
t138 = t153 * t198 - t157 * t154;
t139 = t157 * t153 + t154 * t198;
t122 = t138 * t158 + t139 * t155;
t116 = 0.1e1 / t122;
t200 = t157 * t156;
t142 = atan2(-t200, -t159);
t141 = cos(t142);
t140 = sin(t142);
t186 = t140 * t200;
t126 = -t141 * t159 - t186;
t123 = 0.1e1 / t126;
t149 = 0.1e1 / t159;
t117 = 0.1e1 / t122 ^ 2;
t124 = 0.1e1 / t126 ^ 2;
t218 = -0.2e1 * t156;
t148 = t157 ^ 2;
t145 = t148 * t204 + 0.1e1;
t143 = 0.1e1 / t145;
t217 = t143 - 0.1e1;
t196 = qJD(1) * t160;
t182 = t156 * t196;
t194 = qJD(2) * t159;
t195 = qJD(2) * t157;
t111 = (-(-t157 * t194 - t182) * t149 + t195 * t204) * t143;
t206 = t141 * t156;
t102 = (-t111 * t157 + qJD(2)) * t206 + (-t182 + (t111 - t195) * t159) * t140;
t216 = t102 * t123 * t124;
t121 = t138 * t155 - t139 * t158;
t199 = t157 * t159;
t136 = -t153 * t199 - t154 * t160;
t193 = qJD(2) * t160;
t180 = t156 * t193;
t129 = t136 * qJD(1) - t153 * t180;
t137 = t153 * t160 - t154 * t199;
t130 = t137 * qJD(1) - t154 * t180;
t106 = -t121 * qJD(6) + t129 * t158 + t130 * t155;
t118 = t116 * t117;
t215 = t106 * t118;
t105 = t122 * qJD(6) + t129 * t155 - t130 * t158;
t115 = t121 ^ 2;
t110 = t115 * t117 + 0.1e1;
t211 = t117 * t121;
t214 = 0.1e1 / t110 ^ 2 * (t105 * t211 - t115 * t215);
t213 = t111 * t140;
t212 = t111 * t156;
t170 = -t153 * t158 - t154 * t155;
t201 = t156 * t160;
t134 = t170 * t201;
t210 = t117 * t134;
t209 = t124 * t156;
t202 = t149 * t156;
t168 = qJD(2) * (t149 * t221 + t202);
t172 = t147 * t157 * t196;
t208 = (t148 * t168 + t150 * t172) / t145 ^ 2;
t177 = 0.1e1 + t204;
t128 = t177 * t157 * t143;
t207 = t128 * t157;
t205 = t147 * t149;
t152 = t160 ^ 2;
t203 = t147 * t152;
t114 = t124 * t203 + 0.1e1;
t191 = 0.2e1 * (-t203 * t216 + (t152 * t156 * t194 - t172) * t124) / t114 ^ 2;
t190 = 0.2e1 * t216;
t189 = 0.2e1 * t214;
t188 = 0.2e1 * t118 * t121;
t187 = t124 * t201;
t185 = t143 * t205;
t181 = t156 * t195;
t176 = t156 * t191;
t175 = t208 * t218;
t174 = t208 * t222;
t173 = t157 * t185;
t171 = t177 * t160;
t120 = t136 * t158 + t137 * t155;
t119 = t136 * t155 - t137 * t158;
t169 = -t153 * t155 + t154 * t158;
t133 = t169 * t201;
t132 = -t139 * qJD(1) + t154 * t181;
t131 = -t138 * qJD(1) + t153 * t181;
t112 = 0.1e1 / t114;
t108 = 0.1e1 / t110;
t107 = (t217 * t156 * t140 - t141 * t173) * t160;
t104 = -t140 * t199 + t206 + (t140 * t159 - t141 * t200) * t128;
t103 = -t177 * t174 + (qJD(1) * t171 + t168 * t222) * t143;
t1 = [t160 * t149 * t175 + (qJD(2) * t171 - t197 * t202) * t143, t103, 0, 0, 0, 0; (t123 * t176 + (-t123 * t194 + (qJD(1) * t107 + t102) * t209) * t112) * t157 + (t124 * t176 * t107 + (-((t111 * t173 + t217 * t194 + t175) * t140 + (t174 * t205 - t212 + (t212 + (t218 - t221) * t195) * t143) * t141) * t187 + (-t124 * t194 + t156 * t190) * t107 + (-t123 + ((-t148 + t152) * t141 * t185 + t217 * t186) * t124) * t156 * qJD(1)) * t112) * t160 (t104 * t209 - t123 * t159) * t160 * t191 + ((-t123 * t197 + (-qJD(2) * t104 - t102) * t160 * t124) * t159 + (-t123 * t193 - (-t103 * t141 * t157 + t140 * t195 + t207 * t213 - t213 + (-qJD(2) * t140 - t141 * t196) * t128) * t187 + (t124 * t197 + t160 * t190) * t104 - ((t103 - t196) * t140 + ((0.1e1 - t207) * qJD(2) + (t128 - t157) * t111) * t141) * t124 * t198) * t156) * t112, 0, 0, 0, 0; (-t116 * t119 + t120 * t211) * t189 + ((t120 * qJD(6) + t131 * t155 - t132 * t158) * t116 + t120 * t106 * t188 + (-t119 * t106 - (-t119 * qJD(6) + t131 * t158 + t132 * t155) * t121 - t120 * t105) * t117) * t108 (-t116 * t133 + t121 * t210) * t189 + (-t105 * t210 + (-t117 * t133 + t134 * t188) * t106 + (t169 * t116 - t170 * t211) * t159 * t193 + ((-t219 * t116 - t220 * t211) * t154 + (t220 * t116 - t219 * t211) * t153) * t156) * t108, 0, 0, 0, -0.2e1 * t214 + 0.2e1 * (t105 * t117 * t108 + (-t108 * t215 - t117 * t214) * t121) * t121;];
JaD_rot  = t1;
