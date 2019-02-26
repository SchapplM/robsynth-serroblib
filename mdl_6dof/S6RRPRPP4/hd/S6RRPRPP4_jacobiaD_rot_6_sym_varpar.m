% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPP4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:36:36
% EndTime: 2019-02-26 21:36:37
% DurationCPUTime: 1.02s
% Computational Cost: add. (4821->118), mult. (6168->265), div. (1114->15), fcn. (7752->9), ass. (0->113)
t152 = qJ(4) + pkin(9);
t150 = sin(t152);
t158 = sin(qJ(2));
t161 = cos(qJ(1));
t151 = cos(t152);
t159 = sin(qJ(1));
t210 = t159 * t151;
t139 = t150 * t161 + t158 * t210;
t160 = cos(qJ(2));
t209 = t160 * t151;
t127 = atan2(t139, t209);
t123 = sin(t127);
t124 = cos(t127);
t136 = t139 ^ 2;
t148 = 0.1e1 / t151 ^ 2;
t155 = 0.1e1 / t160 ^ 2;
t215 = t148 * t155;
t128 = t136 * t215 + 0.1e1;
t125 = 0.1e1 / t128;
t147 = 0.1e1 / t151;
t188 = t147 * t155 * t158;
t175 = t139 * t188 + t159;
t114 = t175 * t125;
t207 = t114 - t159;
t234 = t207 * t123 * t160 + t124 * t158;
t211 = t158 * t161;
t138 = t150 * t211 + t210;
t134 = 0.1e1 / t138 ^ 2;
t153 = t160 ^ 2;
t157 = t161 ^ 2;
t212 = t153 * t157;
t191 = t134 * t212;
t131 = 0.1e1 + t191;
t203 = qJD(2) * t160;
t205 = qJD(1) * t161;
t171 = -t153 * t159 * t205 - t157 * t158 * t203;
t179 = qJD(1) * t158 + qJD(4);
t180 = qJD(4) * t158 + qJD(1);
t202 = qJD(2) * t161;
t184 = t160 * t202;
t213 = t151 * t161;
t122 = t180 * t213 + (-t159 * t179 + t184) * t150;
t133 = 0.1e1 / t138;
t224 = t122 * t133 * t134;
t178 = t212 * t224;
t233 = (t134 * t171 - t178) / t131 ^ 2;
t154 = 0.1e1 / t160;
t214 = t150 * t159;
t140 = -t158 * t214 + t213;
t216 = t148 * t150;
t173 = t139 * t216 + t140 * t147;
t232 = t154 * t173;
t208 = t160 * t161;
t222 = t123 * t139;
t118 = t124 * t209 + t222;
t115 = 0.1e1 / t118;
t116 = 0.1e1 / t118 ^ 2;
t187 = t151 * t211;
t137 = -t187 + t214;
t230 = 0.2e1 * t137;
t132 = t137 ^ 2;
t113 = t116 * t132 + 0.1e1;
t121 = qJD(1) * t139 + qJD(4) * t138 - t151 * t184;
t225 = t121 * t116;
t199 = qJD(4) * t160;
t204 = qJD(2) * t158;
t172 = -t150 * t199 - t151 * t204;
t189 = t139 * t215;
t176 = t180 * t159;
t185 = t159 * t203;
t200 = qJD(4) * t151;
t119 = -qJD(1) * t187 + t150 * t176 - t151 * t185 - t161 * t200;
t217 = t147 * t154;
t192 = t119 * t217;
t107 = (-t172 * t189 - t192) * t125;
t170 = -t107 * t139 - t172;
t103 = (-t107 * t209 - t119) * t123 - t170 * t124;
t117 = t115 * t116;
t228 = t103 * t117;
t229 = (-t132 * t228 + t137 * t225) / t113 ^ 2;
t149 = t147 * t148;
t156 = t154 / t153;
t201 = qJD(4) * t150;
t183 = t155 * t201;
t227 = (-t119 * t189 + (t148 * t156 * t204 + t149 * t183) * t136) / t128 ^ 2;
t226 = t116 * t137;
t223 = t123 * t137;
t220 = t124 * t137;
t219 = t124 * t139;
t206 = qJD(1) * t159;
t198 = 0.2e1 * t229;
t197 = -0.2e1 * t227;
t196 = 0.2e1 * t233;
t195 = t117 * t230;
t194 = t115 * t229;
t193 = t116 * t223;
t190 = t139 * t217;
t186 = t155 * t204;
t182 = t116 * t198;
t181 = t208 * t230;
t177 = 0.2e1 * t217 * t227;
t174 = -t134 * t140 * t161 - t133 * t159;
t169 = t121 * t217 - (-t148 * t154 * t201 - t147 * t186) * t137;
t129 = 0.1e1 / t131;
t120 = t151 * t176 + (t161 * t179 + t185) * t150;
t111 = 0.1e1 / t113;
t110 = t125 * t232;
t106 = (-t123 + (-t124 * t190 + t123) * t125) * t137;
t105 = t114 * t219 - t234 * t151;
t104 = -t124 * t150 * t160 + t123 * t140 + (-t123 * t209 + t219) * t110;
t102 = t175 * t197 + (-t119 * t188 + t205 + (t148 * t158 * t183 + (0.2e1 * t156 * t158 ^ 2 + t154) * t147 * qJD(2)) * t139) * t125;
t100 = t197 * t232 + (t173 * t186 + (-t119 * t216 - t120 * t147 + (t140 * t216 + (0.2e1 * t149 * t150 ^ 2 + t147) * t139) * qJD(4)) * t154) * t125;
t1 = [-t125 * t169 + t137 * t177, t102, 0, t100, 0, 0; -0.2e1 * t139 * t194 + (-t119 * t115 + (-t103 * t139 - t106 * t121) * t116) * t111 + (t106 * t182 + (0.2e1 * t106 * t228 + (-t121 * t125 + t121 - (t107 * t125 * t190 + t197) * t137) * t116 * t123 + (-(t139 * t177 - t107) * t226 + (-(t107 + t192) * t137 + t169 * t139) * t116 * t125) * t124) * t111) * t137, t105 * t137 * t182 + (-(t102 * t219 + (-t107 * t222 - t119 * t124) * t114) * t226 + (t103 * t195 - t225) * t105 + (t115 * t208 - t234 * t226) * t201) * t111 + (0.2e1 * t194 * t208 + ((t115 * t202 - (qJD(2) * t207 + t107) * t193) * t158 + (t115 * t206 + (t161 * t103 - (-t102 + t205) * t223 - (-t107 * t207 - qJD(2)) * t220) * t116) * t160) * t111) * t151, 0 (t104 * t226 - t115 * t138) * t198 + (-t104 * t225 + t122 * t115 + (t104 * t195 - t138 * t116) * t103 - (t150 * t204 - t151 * t199 + t100 * t139 - t110 * t119 + (-t110 * t209 + t140) * t107) * t116 * t220 - (-t120 + (-t100 * t151 + t107 * t150) * t160 + t170 * t110) * t193) * t111, 0, 0; t174 * t160 * t196 + (t174 * t204 + ((qJD(1) * t133 - 0.2e1 * t140 * t224) * t161 + (-t120 * t161 + (-qJD(1) * t140 - t122) * t159) * t134) * t160) * t129 (-t133 * t211 - t150 * t191) * t196 + (-0.2e1 * t150 * t178 + (-t158 * t206 + t184) * t133 + (-t122 * t211 + 0.2e1 * t150 * t171 + t200 * t212) * t134) * t129, 0, t134 * t181 * t233 + (t181 * t224 + (-t121 * t208 + (t158 * t202 + t160 * t206) * t137) * t134) * t129, 0, 0;];
JaD_rot  = t1;
