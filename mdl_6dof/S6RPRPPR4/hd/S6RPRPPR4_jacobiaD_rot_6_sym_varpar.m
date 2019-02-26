% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:54
% EndTime: 2019-02-26 20:40:55
% DurationCPUTime: 0.84s
% Computational Cost: add. (2357->105), mult. (3345->232), div. (468->12), fcn. (4023->11), ass. (0->107)
t165 = pkin(9) + qJ(3);
t163 = sin(t165);
t159 = t163 ^ 2;
t164 = cos(t165);
t161 = 0.1e1 / t164 ^ 2;
t220 = t159 * t161;
t171 = sin(qJ(1));
t166 = t171 ^ 2;
t157 = t166 * t220 + 0.1e1;
t160 = 0.1e1 / t164;
t217 = t160 * t163;
t239 = t163 * t220;
t181 = qJD(3) * (t160 * t239 + t217);
t173 = cos(qJ(1));
t209 = qJD(1) * t173;
t218 = t159 * t171;
t187 = t209 * t218;
t225 = (t161 * t187 + t166 * t181) / t157 ^ 2;
t240 = -0.2e1 * t225;
t191 = 0.1e1 + t220;
t238 = t171 * t191;
t154 = 0.1e1 / t157;
t198 = t163 * t209;
t207 = qJD(3) * t171;
t123 = ((t164 * t207 + t198) * t160 + t207 * t220) * t154;
t237 = t123 - t207;
t170 = sin(qJ(6));
t172 = cos(qJ(6));
t205 = qJD(6) * t173;
t210 = qJD(1) * t171;
t236 = t170 * t210 - t172 * t205;
t235 = t170 * t205 + t172 * t210;
t213 = t171 * t163;
t156 = atan2(t213, t164);
t153 = cos(t156);
t152 = sin(t156);
t200 = t152 * t213;
t130 = t153 * t164 + t200;
t127 = 0.1e1 / t130;
t169 = cos(pkin(10));
t211 = t171 * t169;
t168 = sin(pkin(10));
t215 = t168 * t173;
t150 = t164 * t215 - t211;
t212 = t171 * t168;
t214 = t169 * t173;
t151 = t164 * t214 + t212;
t140 = t150 * t170 + t151 * t172;
t134 = 0.1e1 / t140;
t128 = 0.1e1 / t130 ^ 2;
t135 = 0.1e1 / t140 ^ 2;
t234 = 0.2e1 * t163;
t233 = t154 - 0.1e1;
t167 = t173 ^ 2;
t219 = t159 * t167;
t126 = t128 * t219 + 0.1e1;
t208 = qJD(3) * t164;
t221 = t153 * t163;
t114 = (t123 * t171 - qJD(3)) * t221 + (-t237 * t164 + t198) * t152;
t231 = t114 * t127 * t128;
t232 = (-t219 * t231 + (t163 * t167 * t208 - t187) * t128) / t126 ^ 2;
t148 = -t164 * t212 - t214;
t206 = qJD(3) * t173;
t194 = t163 * t206;
t141 = qJD(1) * t148 - t168 * t194;
t149 = -t164 * t211 + t215;
t142 = qJD(1) * t149 - t169 * t194;
t184 = t150 * t172 - t151 * t170;
t118 = qJD(6) * t184 + t141 * t170 + t142 * t172;
t136 = t134 * t135;
t230 = t118 * t136;
t117 = qJD(6) * t140 - t141 * t172 + t142 * t170;
t133 = t184 ^ 2;
t122 = t133 * t135 + 0.1e1;
t224 = t135 * t184;
t229 = 0.1e1 / t122 ^ 2 * (-t117 * t224 - t133 * t230);
t228 = t123 * t163;
t227 = t128 * t163;
t226 = t128 * t173;
t182 = -t168 * t170 - t169 * t172;
t216 = t163 * t173;
t146 = t182 * t216;
t223 = t135 * t146;
t222 = t152 * t171;
t204 = -0.2e1 * t231;
t203 = 0.2e1 * t229;
t202 = -0.2e1 * t136 * t184;
t201 = t128 * t216;
t199 = t154 * t159 * t160;
t195 = t163 * t207;
t190 = -0.2e1 * t163 * t232;
t189 = t160 * t240;
t188 = t171 * t199;
t186 = t191 * t173;
t185 = t148 * t172 - t149 * t170;
t138 = t148 * t170 + t149 * t172;
t183 = t168 * t172 - t169 * t170;
t145 = t183 * t216;
t144 = -qJD(1) * t151 + t169 * t195;
t143 = -qJD(1) * t150 + t168 * t195;
t132 = t154 * t238;
t124 = 0.1e1 / t126;
t120 = 0.1e1 / t122;
t119 = (-t152 * t163 * t233 + t153 * t188) * t173;
t116 = t164 * t222 - t221 + (-t152 * t164 + t153 * t213) * t132;
t115 = t238 * t240 + (qJD(1) * t186 + 0.2e1 * t171 * t181) * t154;
t1 = [t189 * t216 + (qJD(3) * t186 - t210 * t217) * t154, 0, t115, 0, 0, 0; (t127 * t190 + (t127 * t208 + (-qJD(1) * t119 - t114) * t227) * t124) * t171 + (t128 * t190 * t119 + (((-t123 * t188 - t208 * t233 + t225 * t234) * t152 + (t189 * t218 + t228 + (-t228 + (t234 + t239) * t207) * t154) * t153) * t201 + (t128 * t208 + t163 * t204) * t119 + (t127 + ((-t166 + t167) * t153 * t199 + t233 * t200) * t128) * t163 * qJD(1)) * t124) * t173, 0, 0.2e1 * (-t116 * t227 + t127 * t164) * t173 * t232 + ((t127 * t210 + (qJD(3) * t116 + t114) * t226) * t164 + (t127 * t206 + (t115 * t153 * t171 + t237 * t152 + (qJD(3) * t152 - t123 * t222 + t153 * t209) * t132) * t201 + (-t128 * t210 + t173 * t204) * t116 + ((-t115 + t209) * t152 + ((t132 * t171 - 0.1e1) * qJD(3) + (-t132 + t171) * t123) * t153) * t164 * t226) * t163) * t124, 0, 0, 0; (t134 * t185 - t138 * t224) * t203 + ((qJD(6) * t138 - t143 * t172 + t144 * t170) * t134 + t138 * t118 * t202 + (t185 * t118 + (qJD(6) * t185 + t143 * t170 + t144 * t172) * t184 - t138 * t117) * t135) * t120, 0 (-t134 * t145 - t184 * t223) * t203 + (-t117 * t223 + (-t135 * t145 + t146 * t202) * t118 + (t134 * t183 + t182 * t224) * t164 * t206 + ((t236 * t134 + t235 * t224) * t169 + (-t235 * t134 + t236 * t224) * t168) * t163) * t120, 0, 0, -0.2e1 * t229 - 0.2e1 * (t117 * t135 * t120 - (-t120 * t230 - t135 * t229) * t184) * t184;];
JaD_rot  = t1;
