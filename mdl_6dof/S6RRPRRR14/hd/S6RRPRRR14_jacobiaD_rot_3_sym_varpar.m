% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRR14
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR14_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_rot_3_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:55:55
% EndTime: 2019-02-26 22:55:56
% DurationCPUTime: 0.74s
% Computational Cost: add. (2253->89), mult. (6986->213), div. (424->12), fcn. (8980->13), ass. (0->101)
t171 = sin(pkin(6));
t173 = cos(pkin(7));
t174 = cos(pkin(6));
t170 = sin(pkin(7));
t177 = cos(qJ(2));
t217 = t170 * t177;
t159 = -t171 * t217 + t174 * t173;
t156 = 0.1e1 / t159;
t175 = sin(qJ(2));
t178 = cos(qJ(1));
t210 = t178 * t175;
t176 = sin(qJ(1));
t211 = t176 * t177;
t191 = t174 * t210 + t211;
t216 = t171 * t175;
t157 = 0.1e1 / t159 ^ 2;
t209 = t178 * t177;
t212 = t176 * t175;
t160 = -t174 * t209 + t212;
t214 = t171 * t178;
t195 = -t160 * t170 + t173 * t214;
t220 = t195 * t157;
t231 = t170 * (t156 * t191 + t216 * t220);
t143 = atan2(t195, t159);
t138 = sin(t143);
t139 = cos(t143);
t122 = t138 * t195 + t139 * t159;
t119 = 0.1e1 / t122;
t169 = sin(pkin(14));
t172 = cos(pkin(14));
t190 = t174 * t212 - t209;
t192 = t174 * t211 + t210;
t215 = t171 * t176;
t193 = t170 * t215 - t173 * t192;
t135 = t193 * t169 - t172 * t190;
t129 = 0.1e1 / t135;
t120 = 0.1e1 / t122 ^ 2;
t130 = 0.1e1 / t135 ^ 2;
t154 = -t170 * t192 - t173 * t215;
t151 = t154 ^ 2;
t117 = t151 * t120 + 0.1e1;
t146 = t160 * qJD(1) + t190 * qJD(2);
t208 = qJD(1) * t171;
t201 = t178 * t208;
t136 = t146 * t170 - t173 * t201;
t224 = t136 * t120;
t150 = t195 ^ 2;
t142 = t150 * t157 + 0.1e1;
t140 = 0.1e1 / t142;
t148 = t192 * qJD(1) + t191 * qJD(2);
t202 = t176 * t208;
t137 = -t148 * t170 - t173 * t202;
t207 = qJD(2) * t171;
t218 = t170 * t175;
t198 = t207 * t218;
t197 = t157 * t198;
t185 = t137 * t156 - t195 * t197;
t113 = t185 * t140;
t196 = -t138 * t159 + t139 * t195;
t109 = t196 * t113 + t138 * t137 + t139 * t198;
t229 = t109 * t119 * t120;
t230 = (-t151 * t229 + t154 * t224) / t117 ^ 2;
t158 = t156 * t157;
t228 = (-t150 * t158 * t198 + t137 * t220) / t142 ^ 2;
t227 = t120 * t154;
t147 = t191 * qJD(1) + t192 * qJD(2);
t188 = t146 * t173 + t170 * t201;
t127 = -t147 * t172 + t188 * t169;
t226 = t127 * t129 * t130;
t134 = -t169 * t190 - t193 * t172;
t225 = t130 * t134;
t223 = t138 * t154;
t222 = t139 * t154;
t221 = t195 * t156;
t219 = t169 * t173;
t213 = t172 * t173;
t206 = -0.2e1 * t230;
t205 = -0.2e1 * t229;
t128 = t134 ^ 2;
t125 = t128 * t130 + 0.1e1;
t126 = -t147 * t169 - t188 * t172;
t204 = 0.2e1 * (t126 * t225 - t128 * t226) / t125 ^ 2;
t203 = 0.2e1 * t228;
t200 = -0.2e1 * t156 * t228;
t199 = 0.2e1 * t134 * t226;
t194 = t160 * t173 + t170 * t214;
t187 = -t148 * t173 + t170 * t202;
t186 = t138 + (t139 * t221 - t138) * t140;
t168 = t170 ^ 2;
t149 = t190 * qJD(1) + t160 * qJD(2);
t145 = -t172 * t192 + t190 * t219;
t144 = -t169 * t192 - t190 * t213;
t133 = t194 * t169 - t172 * t191;
t132 = -t169 * t191 - t194 * t172;
t123 = 0.1e1 / t125;
t115 = 0.1e1 / t117;
t114 = t140 * t231;
t112 = t186 * t154;
t110 = (-t138 * t191 + t139 * t216) * t170 - t196 * t114;
t108 = t203 * t231 + (t149 * t156 * t170 + (-t137 * t157 * t218 + (t157 * t191 * t168 * t175 + (0.2e1 * t158 * t168 * t171 * t175 ^ 2 - t157 * t217) * t195) * qJD(2)) * t171) * t140;
t1 = [t154 * t200 + (t136 * t156 - t154 * t197) * t140, t108, 0, 0, 0, 0; t195 * t119 * t206 + (t137 * t119 + (-t109 * t195 + t112 * t136) * t120) * t115 + ((t112 * t205 + t186 * t224) * t115 + (t112 * t206 + ((-t113 * t140 * t221 + t203) * t223 + (t195 * t200 + t113 + (-t113 + t185) * t140) * t222) * t115) * t120) * t154, 0.2e1 * (t119 * t170 * t190 - t110 * t227) * t230 + ((t196 * t108 - (-t122 * t113 + t137 * t139) * t114) * t227 + (t154 * t205 + t224) * t110 + (-t147 * t119 + (t190 * t109 + (-t113 * t191 + t177 * t207) * t222 + (t149 + (qJD(2) * t114 - t113) * t216) * t223) * t120) * t170) * t115, 0, 0, 0, 0; (-t129 * t132 + t133 * t225) * t204 + ((t149 * t169 + t187 * t172) * t129 + t133 * t199 + (-t132 * t127 - (t149 * t172 - t187 * t169) * t134 - t133 * t126) * t130) * t123 (-t129 * t144 + t145 * t225) * t204 + ((t146 * t169 - t147 * t213) * t129 + t145 * t199 + (-t144 * t127 - (t146 * t172 + t147 * t219) * t134 - t145 * t126) * t130) * t123, 0, 0, 0, 0;];
JaD_rot  = t1;
