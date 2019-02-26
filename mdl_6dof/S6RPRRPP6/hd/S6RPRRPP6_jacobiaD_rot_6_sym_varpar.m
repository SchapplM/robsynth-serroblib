% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:59:25
% EndTime: 2019-02-26 20:59:26
% DurationCPUTime: 0.96s
% Computational Cost: add. (4821->119), mult. (6168->264), div. (1114->15), fcn. (7752->9), ass. (0->112)
t159 = sin(qJ(1));
t160 = cos(qJ(3));
t209 = t159 * t160;
t152 = qJ(4) + pkin(9);
t150 = sin(t152);
t158 = sin(qJ(3));
t230 = cos(qJ(1));
t191 = t230 * t158;
t151 = cos(t152);
t210 = t159 * t151;
t139 = t150 * t191 + t210;
t208 = t160 * t150;
t127 = atan2(t139, t208);
t124 = cos(t127);
t123 = sin(t127);
t221 = t123 * t139;
t118 = t124 * t208 + t221;
t115 = 0.1e1 / t118;
t138 = t150 * t230 + t158 * t210;
t133 = 0.1e1 / t138;
t147 = 0.1e1 / t150;
t155 = 0.1e1 / t160;
t116 = 0.1e1 / t118 ^ 2;
t134 = 0.1e1 / t138 ^ 2;
t148 = 0.1e1 / t150 ^ 2;
t211 = t159 * t150;
t137 = -t151 * t230 + t158 * t211;
t132 = t137 ^ 2;
t113 = t116 * t132 + 0.1e1;
t186 = qJD(1) * t230;
t204 = qJD(3) * t160;
t172 = -t158 * t186 - t159 * t204;
t169 = qJD(4) * t230 - t172;
t181 = qJD(4) * t158 + qJD(1);
t121 = t150 * t169 + t181 * t210;
t224 = t121 * t116;
t136 = t139 ^ 2;
t156 = 0.1e1 / t160 ^ 2;
t214 = t148 * t156;
t128 = t136 * t214 + 0.1e1;
t125 = 0.1e1 / t128;
t202 = qJD(4) * t160;
t206 = qJD(3) * t158;
t173 = -t150 * t206 + t151 * t202;
t193 = t139 * t214;
t190 = t230 * t160;
t178 = qJD(3) * t190;
t179 = t151 * t191;
t119 = -qJD(4) * t179 - t150 * t178 - t151 * t186 + (qJD(1) * t158 + qJD(4)) * t211;
t216 = t147 * t155;
t196 = t119 * t216;
t107 = (-t173 * t193 - t196) * t125;
t171 = -t107 * t139 - t173;
t103 = (-t107 * t208 - t119) * t123 - t171 * t124;
t117 = t115 * t116;
t228 = t103 * t117;
t229 = (-t132 * t228 + t137 * t224) / t113 ^ 2;
t149 = t147 * t148;
t154 = t160 ^ 2;
t157 = t155 / t154;
t203 = qJD(4) * t151;
t188 = t156 * t203;
t227 = (-t119 * t193 + (t148 * t157 * t206 - t149 * t188) * t136) / t128 ^ 2;
t153 = t159 ^ 2;
t213 = t153 * t154;
t195 = t134 * t213;
t131 = 0.1e1 + t195;
t170 = -t153 * t158 * t204 + t154 * t159 * t186;
t122 = t151 * t169 - t181 * t211;
t223 = t122 * t133 * t134;
t180 = t213 * t223;
t226 = (t134 * t170 - t180) / t131 ^ 2;
t225 = t116 * t137;
t222 = t123 * t137;
t220 = t123 * t160;
t219 = t124 * t137;
t218 = t124 * t139;
t217 = t124 * t158;
t215 = t148 * t151;
t212 = t158 * t159;
t207 = qJD(1) * t159;
t205 = qJD(3) * t159;
t201 = 0.2e1 * t229;
t200 = -0.2e1 * t227;
t199 = 0.2e1 * t226;
t198 = 0.2e1 * t117 * t137;
t197 = t116 * t222;
t194 = t139 * t216;
t192 = t147 * t156 * t158;
t189 = t156 * t206;
t175 = t139 * t192 + t230;
t114 = t175 * t125;
t187 = t230 - t114;
t185 = -0.2e1 * t115 * t229;
t184 = t116 * t201;
t183 = 0.2e1 * t155 * t227;
t182 = -0.2e1 * t137 * t209;
t177 = t147 * t183;
t140 = t179 - t211;
t176 = t139 * t215 - t140 * t147;
t174 = t134 * t140 * t159 - t133 * t230;
t168 = t121 * t216 - (t148 * t155 * t203 - t147 * t189) * t137;
t129 = 0.1e1 / t131;
t120 = qJD(1) * t138 + qJD(4) * t139 - t151 * t178;
t111 = 0.1e1 / t113;
t110 = t176 * t155 * t125;
t106 = (-t123 + (-t124 * t194 + t123) * t125) * t137;
t105 = t114 * t218 + (t187 * t220 - t217) * t150;
t104 = t124 * t151 * t160 + t123 * t140 - (-t123 * t208 + t218) * t110;
t102 = t175 * t200 + (-t119 * t192 - t207 + (-t148 * t158 * t188 + (0.2e1 * t157 * t158 ^ 2 + t155) * t147 * qJD(3)) * t139) * t125;
t100 = t176 * t183 + (-t176 * t189 + (t119 * t215 - t120 * t147 + (-t140 * t215 + (0.2e1 * t149 * t151 ^ 2 + t147) * t139) * qJD(4)) * t155) * t125;
t1 = [-t125 * t168 + t137 * t177, 0, t102, t100, 0, 0; t139 * t185 + (-t119 * t115 + (-t103 * t139 - t106 * t121) * t116) * t111 + (t106 * t184 + (0.2e1 * t106 * t228 + (-t121 * t125 + t121 - (t107 * t125 * t194 + t200) * t137) * t116 * t123 + (-(t139 * t177 - t107) * t225 + (-(t107 + t196) * t137 + t168 * t139) * t116 * t125) * t124) * t111) * t137, 0, t105 * t137 * t184 + (-(t102 * t218 + (-t107 * t221 - t119 * t124) * t114) * t225 + (t103 * t198 - t224) * t105 + (t115 * t209 - (-t114 * t220 + t123 * t190 - t217) * t225) * t203) * t111 + (t185 * t209 + ((-t115 * t205 - (-qJD(3) * t187 + t107) * t197) * t158 + (t115 * t186 + (-t159 * t103 - (-t102 - t207) * t222 - (t107 * t187 - qJD(3)) * t219) * t116) * t160) * t111) * t150 (t104 * t225 - t115 * t138) * t201 + (-t104 * t224 + t122 * t115 + (t104 * t198 - t138 * t116) * t103 - (-t151 * t206 - t150 * t202 + t100 * t139 + t110 * t119 + (t110 * t208 + t140) * t107) * t116 * t219 - (-t120 + (-t100 * t150 - t107 * t151) * t160 - t171 * t110) * t197) * t111, 0, 0; t174 * t160 * t199 + (t174 * t206 + ((-qJD(1) * t133 + 0.2e1 * t140 * t223) * t159 + (t120 * t159 - t122 * t230 - t140 * t186) * t134) * t160) * t129, 0 (t133 * t212 + t151 * t195) * t199 + (0.2e1 * t151 * t180 + t172 * t133 + (qJD(4) * t150 * t213 + t122 * t212 - 0.2e1 * t151 * t170) * t134) * t129, t134 * t182 * t226 + (t182 * t223 + (t121 * t209 + (-t158 * t205 + t160 * t186) * t137) * t134) * t129, 0, 0;];
JaD_rot  = t1;
