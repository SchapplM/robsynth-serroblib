% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR10_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:10
% EndTime: 2019-02-26 21:43:11
% DurationCPUTime: 0.69s
% Computational Cost: add. (1788->92), mult. (4303->201), div. (668->14), fcn. (5516->11), ass. (0->91)
t178 = sin(qJ(2));
t179 = sin(qJ(1));
t230 = cos(pkin(6));
t201 = t179 * t230;
t199 = t178 * t201;
t180 = cos(qJ(2));
t181 = cos(qJ(1));
t215 = t181 * t180;
t162 = -t199 + t215;
t173 = pkin(11) + qJ(4);
t169 = sin(t173);
t170 = cos(t173);
t177 = sin(pkin(6));
t219 = t177 * t179;
t191 = -t162 * t169 + t170 * t219;
t232 = t191 * qJD(4);
t200 = t181 * t230;
t198 = t180 * t200;
t216 = t179 * t178;
t158 = -t198 + t216;
t218 = t177 * t180;
t152 = atan2(-t158, -t218);
t150 = sin(t152);
t151 = cos(t152);
t156 = t158 ^ 2;
t172 = 0.1e1 / t177 ^ 2;
t175 = 0.1e1 / t180 ^ 2;
t155 = t156 * t172 * t175 + 0.1e1;
t153 = 0.1e1 / t155;
t171 = 0.1e1 / t177;
t174 = 0.1e1 / t180;
t205 = t158 * t171 * t174;
t231 = (t151 * t205 - t150) * t153 + t150;
t134 = -t150 * t158 - t151 * t218;
t131 = 0.1e1 / t134;
t149 = t162 * t170 + t169 * t219;
t143 = 0.1e1 / t149;
t132 = 0.1e1 / t134 ^ 2;
t144 = 0.1e1 / t149 ^ 2;
t188 = -t178 * t200 - t179 * t180;
t189 = -t181 * t178 - t180 * t201;
t140 = -t189 * qJD(1) - t188 * qJD(2);
t213 = qJD(2) * t178;
t202 = t175 * t213;
t190 = t140 * t174 + t158 * t202;
t221 = t153 * t171;
t123 = t190 * t221;
t194 = t150 * t218 - t151 * t158;
t206 = t151 * t177 * t178;
t119 = qJD(2) * t206 + t194 * t123 - t150 * t140;
t229 = t119 * t131 * t132;
t139 = t188 * qJD(1) + t189 * qJD(2);
t214 = qJD(1) * t177;
t203 = t181 * t214;
t128 = t149 * qJD(4) + t139 * t169 - t170 * t203;
t142 = t191 ^ 2;
t137 = t142 * t144 + 0.1e1;
t224 = t144 * t191;
t129 = t139 * t170 + t169 * t203 + t232;
t227 = t129 * t143 * t144;
t228 = (-t128 * t224 - t142 * t227) / t137 ^ 2;
t176 = t174 * t175;
t226 = (t140 * t158 * t175 + t156 * t176 * t213) * t172 / t155 ^ 2;
t225 = t132 * t189;
t223 = t150 * t189;
t222 = t151 * t189;
t220 = t175 * t178;
t217 = t177 * t181;
t212 = qJD(2) * t180;
t157 = t189 ^ 2;
t127 = t157 * t132 + 0.1e1;
t197 = qJD(2) * t230 + qJD(1);
t138 = -qJD(1) * t198 - t181 * t212 + t197 * t216;
t211 = 0.2e1 * (t138 * t225 - t157 * t229) / t127 ^ 2;
t210 = 0.2e1 * t229;
t209 = 0.2e1 * t228;
t208 = -0.2e1 * t226;
t207 = t191 * t227;
t204 = t179 * t214;
t195 = t169 * t143 + t170 * t224;
t193 = t158 * t220 - t174 * t188;
t192 = -t169 * t188 + t170 * t217;
t147 = t169 * t217 + t170 * t188;
t141 = -qJD(1) * t199 - t179 * t213 + t197 * t215;
t135 = 0.1e1 / t137;
t125 = 0.1e1 / t127;
t124 = t193 * t221;
t122 = t231 * t189;
t120 = t194 * t124 + t150 * t188 + t206;
t118 = (t193 * t208 + (t140 * t220 + t141 * t174 + (-t188 * t220 + (0.2e1 * t176 * t178 ^ 2 + t174) * t158) * qJD(2)) * t153) * t171;
t1 = [(-t189 * t174 * t208 + (-t138 * t174 - t189 * t202) * t153) * t171, t118, 0, 0, 0, 0; t158 * t131 * t211 + (-t140 * t131 + (t119 * t158 + t122 * t138) * t132) * t125 - (t122 * t210 * t125 + (t122 * t211 + ((t123 * t153 * t205 + t208) * t223 + (0.2e1 * t205 * t226 - t123 + (-t190 * t171 + t123) * t153) * t222 - t231 * t138) * t125) * t132) * t189 (-t120 * t225 - t131 * t162) * t211 + (-t120 * t189 * t210 + t139 * t131 + (-t162 * t119 + t120 * t138 + (t177 * t212 - t118 * t158 - t124 * t140 + (t124 * t218 + t188) * t123) * t222 + (t123 * t124 * t158 - t141 + (t118 * t180 + (-qJD(2) * t124 - t123) * t178) * t177) * t223) * t132) * t125, 0, 0, 0, 0; (t143 * t192 - t147 * t224) * t209 + ((t147 * qJD(4) - t141 * t169 + t170 * t204) * t143 - 0.2e1 * t147 * t207 + (t192 * t129 + (t192 * qJD(4) - t141 * t170 - t169 * t204) * t191 - t147 * t128) * t144) * t135, -t195 * t189 * t209 + (t195 * t138 - ((-qJD(4) * t143 + 0.2e1 * t207) * t170 + (t128 * t170 + (t129 + t232) * t169) * t144) * t189) * t135, 0, -0.2e1 * t228 - 0.2e1 * (t128 * t144 * t135 - (-t135 * t227 - t144 * t228) * t191) * t191, 0, 0;];
JaD_rot  = t1;
