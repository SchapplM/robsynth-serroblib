% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR11_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:48
% EndTime: 2019-02-26 21:59:49
% DurationCPUTime: 0.72s
% Computational Cost: add. (2159->94), mult. (2949->208), div. (516->12), fcn. (3430->9), ass. (0->96)
t156 = sin(qJ(2));
t150 = 0.1e1 / t156 ^ 2;
t158 = cos(qJ(2));
t154 = t158 ^ 2;
t204 = t150 * t154;
t222 = t158 * t204;
t157 = sin(qJ(1));
t179 = 0.1e1 + t204;
t221 = t157 * t179;
t147 = qJD(4) + qJD(5) + qJD(6);
t175 = qJD(1) * t156 + t147;
t159 = cos(qJ(1));
t191 = qJD(2) * t159;
t220 = t175 * t157 - t158 * t191;
t192 = qJD(2) * t158;
t219 = t157 * t192 + t175 * t159;
t198 = t157 * t158;
t141 = atan2(-t198, t156);
t139 = cos(t141);
t138 = sin(t141);
t184 = t138 * t198;
t127 = t139 * t156 - t184;
t124 = 0.1e1 / t127;
t148 = qJ(4) + qJ(5) + qJ(6);
t145 = sin(t148);
t146 = cos(t148);
t199 = t157 * t146;
t200 = t156 * t159;
t135 = t145 * t200 + t199;
t131 = 0.1e1 / t135;
t149 = 0.1e1 / t156;
t125 = 0.1e1 / t127 ^ 2;
t132 = 0.1e1 / t135 ^ 2;
t152 = t157 ^ 2;
t144 = t152 * t204 + 0.1e1;
t142 = 0.1e1 / t144;
t218 = t142 - 0.1e1;
t176 = t147 * t156 + qJD(1);
t170 = t176 * t159;
t115 = t145 * t170 + t146 * t220;
t134 = t145 * t157 - t146 * t200;
t130 = t134 ^ 2;
t120 = t130 * t132 + 0.1e1;
t207 = t132 * t134;
t116 = -t145 * t220 + t146 * t170;
t215 = t116 * t131 * t132;
t217 = (t115 * t207 - t130 * t215) / t120 ^ 2;
t195 = qJD(1) * t159;
t182 = t158 * t195;
t193 = qJD(2) * t157;
t117 = ((t156 * t193 - t182) * t149 + t193 * t204) * t142;
t205 = t139 * t158;
t111 = (-t117 * t157 + qJD(2)) * t205 + (-t182 + (-t117 + t193) * t156) * t138;
t216 = t111 * t124 * t125;
t214 = t117 * t138;
t213 = t117 * t158;
t212 = t125 * t158;
t211 = t125 * t159;
t168 = qJD(2) * (-t158 - t222) * t149;
t202 = t154 * t157;
t173 = t195 * t202;
t210 = (t150 * t173 + t152 * t168) / t144 ^ 2;
t129 = t142 * t221;
t209 = t129 * t157;
t208 = t131 * t146;
t206 = t134 * t145;
t155 = t159 ^ 2;
t203 = t154 * t155;
t201 = t156 * t157;
t197 = qJD(1) * t157;
t196 = qJD(1) * t158;
t194 = qJD(2) * t156;
t123 = t125 * t203 + 0.1e1;
t190 = 0.2e1 * (-t203 * t216 + (-t155 * t156 * t192 - t173) * t125) / t123 ^ 2;
t189 = 0.2e1 * t217;
t188 = 0.2e1 * t216;
t187 = -0.2e1 * t210;
t186 = t158 * t211;
t185 = t158 * t210;
t183 = t149 * t202;
t178 = t158 * t190;
t177 = 0.2e1 * t134 * t215;
t174 = t142 * t183;
t172 = t179 * t159;
t171 = t176 * t157;
t169 = t132 * t206 + t208;
t167 = t169 * t159;
t137 = -t145 * t201 + t146 * t159;
t136 = t145 * t159 + t156 * t199;
t121 = 0.1e1 / t123;
t118 = 0.1e1 / t120;
t114 = (t218 * t158 * t138 + t139 * t174) * t159;
t113 = t138 * t201 + t205 + (-t138 * t156 - t139 * t198) * t129;
t112 = t187 * t221 + (qJD(1) * t172 + 0.2e1 * t157 * t168) * t142;
t108 = -0.2e1 * t217 + 0.2e1 * (t115 * t118 * t132 + (-t118 * t215 - t132 * t217) * t134) * t134;
t1 = [0.2e1 * t149 * t159 * t185 + (t149 * t157 * t196 + qJD(2) * t172) * t142, t112, 0, 0, 0, 0; (t124 * t178 + (t124 * t194 + (qJD(1) * t114 + t111) * t212) * t121) * t157 + (t125 * t178 * t114 + (-((-t117 * t174 - t218 * t194 - 0.2e1 * t185) * t138 + (t183 * t187 - t213 + (t213 + (-0.2e1 * t158 - t222) * t193) * t142) * t139) * t186 + (t125 * t194 + t158 * t188) * t114 + (-t124 + ((t152 - t155) * t154 * t149 * t142 * t139 + t218 * t184) * t125) * t196) * t121) * t159 (t113 * t212 + t124 * t156) * t159 * t190 + ((t124 * t197 + (qJD(2) * t113 + t111) * t211) * t156 + (-t124 * t191 - (-t112 * t139 * t157 + t138 * t193 + t209 * t214 - t214 + (-qJD(2) * t138 - t139 * t195) * t129) * t186 + (t125 * t197 + t159 * t188) * t113 - ((-t112 + t195) * t138 + ((-0.1e1 + t209) * qJD(2) + (-t129 + t157) * t117) * t139) * t125 * t200) * t158) * t121, 0, 0, 0, 0; (-t131 * t136 + t137 * t207) * t189 + (t137 * t177 - t131 * t145 * t171 + t219 * t208 + (t134 * t146 * t171 - t137 * t115 - t136 * t116 + t206 * t219) * t132) * t118, t158 * t167 * t189 + (t167 * t194 + (t169 * t197 + ((t131 * t147 + t177) * t145 + (-t115 * t145 + (-t134 * t147 + t116) * t146) * t132) * t159) * t158) * t118, 0, t108, t108, t108;];
JaD_rot  = t1;
