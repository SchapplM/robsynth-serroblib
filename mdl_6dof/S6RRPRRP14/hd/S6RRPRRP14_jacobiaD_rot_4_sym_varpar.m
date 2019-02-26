% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRP14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP14_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:53:24
% EndTime: 2019-02-26 21:53:25
% DurationCPUTime: 0.60s
% Computational Cost: add. (1334->90), mult. (4303->199), div. (668->14), fcn. (5516->11), ass. (0->90)
t174 = sin(qJ(2));
t175 = sin(qJ(1));
t177 = cos(qJ(2));
t178 = cos(qJ(1));
t224 = cos(pkin(6));
t194 = t178 * t224;
t157 = t174 * t194 + t175 * t177;
t172 = sin(pkin(6));
t213 = t172 * t174;
t151 = atan2(-t157, t213);
t147 = sin(t151);
t148 = cos(t151);
t154 = t157 ^ 2;
t168 = 0.1e1 / t172 ^ 2;
t170 = 0.1e1 / t174 ^ 2;
t152 = t154 * t168 * t170 + 0.1e1;
t149 = 0.1e1 / t152;
t167 = 0.1e1 / t172;
t169 = 0.1e1 / t174;
t199 = t157 * t167 * t169;
t225 = (t148 * t199 + t147) * t149 - t147;
t131 = -t147 * t157 + t148 * t213;
t128 = 0.1e1 / t131;
t195 = t175 * t224;
t159 = t178 * t174 + t177 * t195;
t173 = sin(qJ(4));
t176 = cos(qJ(4));
t212 = t172 * t175;
t144 = t159 * t173 + t176 * t212;
t140 = 0.1e1 / t144;
t129 = 0.1e1 / t131 ^ 2;
t141 = 0.1e1 / t144 ^ 2;
t190 = qJD(2) * t224 + qJD(1);
t192 = t174 * t195;
t207 = qJD(2) * t174;
t209 = t178 * t177;
t138 = -qJD(1) * t192 - t175 * t207 + t190 * t209;
t206 = qJD(2) * t177;
t196 = t170 * t206;
t185 = -t138 * t169 + t157 * t196;
t215 = t149 * t167;
t120 = t185 * t215;
t187 = -t147 * t213 - t148 * t157;
t200 = t148 * t172 * t177;
t116 = qJD(2) * t200 + t187 * t120 - t147 * t138;
t223 = t116 * t128 * t129;
t191 = t177 * t194;
t210 = t175 * t174;
t135 = -qJD(1) * t191 - t178 * t206 + t190 * t210;
t208 = qJD(1) * t172;
t197 = t178 * t208;
t126 = t144 * qJD(4) + t135 * t176 + t173 * t197;
t143 = -t159 * t176 + t173 * t212;
t139 = t143 ^ 2;
t134 = t139 * t141 + 0.1e1;
t218 = t141 * t143;
t205 = qJD(4) * t143;
t127 = -t135 * t173 + t176 * t197 - t205;
t220 = t127 * t140 * t141;
t222 = (t126 * t218 - t139 * t220) / t134 ^ 2;
t171 = t169 * t170;
t221 = (t138 * t157 * t170 - t154 * t171 * t206) * t168 / t152 ^ 2;
t160 = -t192 + t209;
t219 = t129 * t160;
t217 = t147 * t160;
t216 = t148 * t160;
t214 = t170 * t177;
t211 = t172 * t178;
t155 = t160 ^ 2;
t124 = t155 * t129 + 0.1e1;
t136 = t157 * qJD(1) + t159 * qJD(2);
t204 = 0.2e1 * (-t136 * t219 - t155 * t223) / t124 ^ 2;
t203 = 0.2e1 * t223;
t202 = 0.2e1 * t222;
t201 = -0.2e1 * t221;
t198 = t175 * t208;
t193 = 0.2e1 * t143 * t220;
t188 = t176 * t140 + t173 * t218;
t156 = -t191 + t210;
t186 = t156 * t169 + t157 * t214;
t146 = -t156 * t173 + t176 * t211;
t145 = t156 * t176 + t173 * t211;
t137 = t159 * qJD(1) + t157 * qJD(2);
t132 = 0.1e1 / t134;
t122 = 0.1e1 / t124;
t121 = t186 * t215;
t119 = t225 * t160;
t117 = t187 * t121 + t147 * t156 + t200;
t115 = (t186 * t201 + (t138 * t214 + t137 * t169 + (-t156 * t214 + (-0.2e1 * t171 * t177 ^ 2 - t169) * t157) * qJD(2)) * t149) * t167;
t1 = [(0.2e1 * t160 * t169 * t221 + (t136 * t169 + t160 * t196) * t149) * t167, t115, 0, 0, 0, 0; t157 * t128 * t204 + (-t138 * t128 + (t116 * t157 + t119 * t136) * t129) * t122 + (t119 * t203 * t122 + (t119 * t204 + (-(-t120 * t149 * t199 + t201) * t217 - (t199 * t201 - t120 + (-t185 * t167 + t120) * t149) * t216 + t225 * t136) * t122) * t129) * t160 (t117 * t219 + t128 * t159) * t204 + (t117 * t160 * t203 + t135 * t128 + (t159 * t116 + t117 * t136 - (-t172 * t207 - t115 * t157 - t121 * t138 + (-t121 * t213 + t156) * t120) * t216 - (t120 * t121 * t157 + t137 + (-t115 * t174 + (-qJD(2) * t121 - t120) * t177) * t172) * t217) * t129) * t122, 0, 0, 0, 0; (-t140 * t145 + t146 * t218) * t202 + ((t146 * qJD(4) + t137 * t176 - t173 * t198) * t140 + t146 * t193 + (-t145 * t127 - (-t145 * qJD(4) - t137 * t173 - t176 * t198) * t143 - t146 * t126) * t141) * t132, t188 * t160 * t202 + (t188 * t136 + ((qJD(4) * t140 + t193) * t173 + (-t126 * t173 + (t127 - t205) * t176) * t141) * t160) * t132, 0, -0.2e1 * t222 + 0.2e1 * (t126 * t141 * t132 + (-t132 * t220 - t141 * t222) * t143) * t143, 0, 0;];
JaD_rot  = t1;
