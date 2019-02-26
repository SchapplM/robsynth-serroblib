% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR7_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:32:07
% EndTime: 2019-02-26 21:32:07
% DurationCPUTime: 0.60s
% Computational Cost: add. (1334->90), mult. (4303->199), div. (668->14), fcn. (5516->11), ass. (0->90)
t170 = sin(qJ(2));
t171 = sin(qJ(1));
t173 = cos(qJ(2));
t174 = cos(qJ(1));
t220 = cos(pkin(6));
t190 = t174 * t220;
t153 = t170 * t190 + t171 * t173;
t168 = sin(pkin(6));
t209 = t168 * t170;
t148 = atan2(-t153, t209);
t144 = sin(t148);
t145 = cos(t148);
t150 = t153 ^ 2;
t164 = 0.1e1 / t168 ^ 2;
t166 = 0.1e1 / t170 ^ 2;
t149 = t150 * t164 * t166 + 0.1e1;
t146 = 0.1e1 / t149;
t163 = 0.1e1 / t168;
t165 = 0.1e1 / t170;
t195 = t153 * t163 * t165;
t221 = (t145 * t195 + t144) * t146 - t144;
t128 = -t144 * t153 + t145 * t209;
t125 = 0.1e1 / t128;
t191 = t171 * t220;
t155 = t174 * t170 + t173 * t191;
t169 = sin(qJ(5));
t172 = cos(qJ(5));
t208 = t168 * t171;
t143 = t155 * t172 - t169 * t208;
t137 = 0.1e1 / t143;
t126 = 0.1e1 / t128 ^ 2;
t138 = 0.1e1 / t143 ^ 2;
t186 = qJD(2) * t220 + qJD(1);
t188 = t170 * t191;
t203 = qJD(2) * t170;
t205 = t174 * t173;
t135 = -qJD(1) * t188 - t171 * t203 + t186 * t205;
t202 = qJD(2) * t173;
t192 = t166 * t202;
t181 = -t135 * t165 + t153 * t192;
t211 = t146 * t163;
t117 = t181 * t211;
t183 = -t144 * t209 - t145 * t153;
t196 = t145 * t168 * t173;
t113 = qJD(2) * t196 + t183 * t117 - t144 * t135;
t219 = t113 * t125 * t126;
t187 = t173 * t190;
t206 = t171 * t170;
t132 = -qJD(1) * t187 - t174 * t202 + t186 * t206;
t204 = qJD(1) * t168;
t193 = t174 * t204;
t123 = t143 * qJD(5) - t132 * t169 + t172 * t193;
t142 = t155 * t169 + t172 * t208;
t136 = t142 ^ 2;
t131 = t136 * t138 + 0.1e1;
t214 = t138 * t142;
t201 = qJD(5) * t142;
t124 = -t132 * t172 - t169 * t193 - t201;
t216 = t124 * t137 * t138;
t218 = (t123 * t214 - t136 * t216) / t131 ^ 2;
t167 = t165 * t166;
t217 = (t135 * t153 * t166 - t150 * t167 * t202) * t164 / t149 ^ 2;
t156 = -t188 + t205;
t215 = t126 * t156;
t213 = t144 * t156;
t212 = t145 * t156;
t210 = t166 * t173;
t207 = t168 * t174;
t151 = t156 ^ 2;
t121 = t151 * t126 + 0.1e1;
t133 = t153 * qJD(1) + t155 * qJD(2);
t200 = 0.2e1 * (-t133 * t215 - t151 * t219) / t121 ^ 2;
t199 = 0.2e1 * t219;
t198 = 0.2e1 * t218;
t197 = -0.2e1 * t217;
t194 = t171 * t204;
t189 = 0.2e1 * t142 * t216;
t184 = -t169 * t137 + t172 * t214;
t152 = -t187 + t206;
t182 = t152 * t165 + t153 * t210;
t140 = -t152 * t169 + t172 * t207;
t141 = -t152 * t172 - t169 * t207;
t134 = t155 * qJD(1) + t153 * qJD(2);
t129 = 0.1e1 / t131;
t119 = 0.1e1 / t121;
t118 = t182 * t211;
t116 = t221 * t156;
t114 = t183 * t118 + t144 * t152 + t196;
t112 = (t182 * t197 + (t135 * t210 + t134 * t165 + (-t152 * t210 + (-0.2e1 * t167 * t173 ^ 2 - t165) * t153) * qJD(2)) * t146) * t163;
t1 = [(0.2e1 * t156 * t165 * t217 + (t133 * t165 + t156 * t192) * t146) * t163, t112, 0, 0, 0, 0; t153 * t125 * t200 + (-t135 * t125 + (t113 * t153 + t116 * t133) * t126) * t119 + (t116 * t199 * t119 + (t116 * t200 + (-(-t117 * t146 * t195 + t197) * t213 - (t195 * t197 - t117 + (-t181 * t163 + t117) * t146) * t212 + t221 * t133) * t119) * t126) * t156 (t114 * t215 + t125 * t155) * t200 + (t114 * t156 * t199 + t132 * t125 + (t155 * t113 + t114 * t133 - (-t168 * t203 - t112 * t153 - t118 * t135 + (-t118 * t209 + t152) * t117) * t212 - (t117 * t118 * t153 + t134 + (-t112 * t170 + (-qJD(2) * t118 - t117) * t173) * t168) * t213) * t126) * t119, 0, 0, 0, 0; (-t137 * t140 + t141 * t214) * t198 + ((t141 * qJD(5) - t134 * t169 - t172 * t194) * t137 + t141 * t189 + (-t140 * t124 - (-t140 * qJD(5) - t134 * t172 + t169 * t194) * t142 - t141 * t123) * t138) * t129, t184 * t156 * t198 + (t184 * t133 + ((qJD(5) * t137 + t189) * t172 + (-t123 * t172 + (-t124 + t201) * t169) * t138) * t156) * t129, 0, 0, -0.2e1 * t218 + 0.2e1 * (t123 * t138 * t129 + (-t129 * t216 - t138 * t218) * t142) * t142, 0;];
JaD_rot  = t1;
