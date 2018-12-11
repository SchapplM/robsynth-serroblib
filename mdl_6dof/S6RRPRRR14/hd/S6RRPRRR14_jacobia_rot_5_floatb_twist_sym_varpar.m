% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR14_jacobia_rot_5_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_rot_5_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_rot_5_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:21
% EndTime: 2018-12-10 18:38:22
% DurationCPUTime: 0.88s
% Computational Cost: add. (9934->92), mult. (9652->166), div. (85->9), fcn. (9681->31), ass. (0->103)
t231 = pkin(6) - qJ(2);
t223 = cos(t231);
t193 = pkin(6) + qJ(2);
t227 = cos(t193) / 0.2e1;
t178 = t223 / 0.2e1 + t227;
t203 = sin(qJ(2));
t204 = sin(qJ(1));
t208 = cos(qJ(1));
t165 = -t208 * t178 + t204 * t203;
t221 = sin(t231);
t226 = sin(t193) / 0.2e1;
t175 = t226 - t221 / 0.2e1;
t207 = cos(qJ(2));
t166 = t208 * t175 + t204 * t207;
t192 = pkin(7) - pkin(14);
t180 = cos(t192) / 0.2e1;
t191 = pkin(7) + pkin(14);
t190 = cos(t191);
t171 = t180 + t190 / 0.2e1;
t194 = sin(pkin(14));
t197 = sin(pkin(6));
t179 = sin(t191) / 0.2e1;
t189 = sin(t192);
t228 = t179 + t189 / 0.2e1;
t219 = t197 * t228;
t154 = t165 * t171 + t166 * t194 + t208 * t219;
t170 = t179 - t189 / 0.2e1;
t198 = cos(pkin(14));
t172 = t180 - t190 / 0.2e1;
t232 = t197 * t172;
t155 = t165 * t170 - t166 * t198 + t208 * t232;
t229 = pkin(8) + qJ(4);
t217 = sin(t229) / 0.2e1;
t230 = pkin(8) - qJ(4);
t220 = sin(t230);
t173 = t217 - t220 / 0.2e1;
t218 = cos(t229) / 0.2e1;
t222 = cos(t230);
t176 = t218 - t222 / 0.2e1;
t206 = cos(qJ(4));
t196 = sin(pkin(7));
t240 = cos(pkin(7));
t224 = t197 * t240;
t211 = t165 * t196 - t208 * t224;
t137 = t154 * t173 + t155 * t206 + t211 * t176;
t202 = sin(qJ(4));
t214 = t217 + t220 / 0.2e1;
t215 = t222 / 0.2e1 + t218;
t241 = -t154 * t215 + t155 * t202 + t211 * t214;
t174 = t226 + t221 / 0.2e1;
t177 = t227 - t223 / 0.2e1;
t200 = cos(pkin(6));
t159 = t174 * t171 + t177 * t194 + t200 * t228;
t160 = t174 * t170 + t200 * t172 - t177 * t198;
t164 = -t174 * t196 + t200 * t240;
t144 = -t159 * t215 + t160 * t202 - t164 * t214;
t129 = atan2(t241, t144);
t126 = sin(t129);
t127 = cos(t129);
t124 = t126 * t241 + t127 * t144;
t123 = 0.1e1 / t124 ^ 2;
t168 = -t204 * t178 - t208 * t203;
t169 = t204 * t175 - t208 * t207;
t156 = t168 * t170 - t169 * t198 + t204 * t232;
t209 = -t168 * t171 - t169 * t194 - t204 * t219;
t210 = t168 * t196 - t204 * t224;
t138 = t156 * t202 + t209 * t215 + t210 * t214;
t239 = t123 * t138;
t238 = t127 * t241;
t139 = t156 * t206 - t209 * t173 + t210 * t176;
t195 = sin(pkin(8));
t199 = cos(pkin(8));
t148 = t209 * t195 - t210 * t199;
t201 = sin(qJ(5));
t205 = cos(qJ(5));
t133 = t139 * t205 + t148 * t201;
t131 = 0.1e1 / t133 ^ 2;
t132 = t139 * t201 - t148 * t205;
t237 = t131 * t132;
t143 = 0.1e1 / t144 ^ 2;
t236 = t241 * t143;
t235 = t138 ^ 2 * t123;
t233 = t169 * t196;
t225 = t132 ^ 2 * t131 + 0.1e1;
t216 = -t126 * t144 + t238;
t212 = t196 * t214;
t158 = t168 * t198 + t169 * t170;
t157 = -t168 * t194 + t169 * t171;
t149 = -t157 * t195 - t199 * t233;
t147 = -t154 * t195 - t199 * t211;
t146 = (t177 * t170 + t174 * t198) * t202 - (t177 * t171 - t174 * t194) * t215 + t177 * t212;
t145 = t159 * t173 + t160 * t206 - t164 * t176;
t142 = 0.1e1 / t144;
t141 = t157 * t173 + t158 * t206 + t176 * t233;
t140 = (-t165 * t198 - t166 * t170) * t202 - (t165 * t194 - t166 * t171) * t215 - t166 * t212;
t130 = 0.1e1 / t133;
t128 = 0.1e1 / (t143 * t241 ^ 2 + 0.1e1);
t125 = 0.1e1 / t225;
t122 = 0.1e1 / t124;
t121 = 0.1e1 / (0.1e1 + t235);
t120 = (-t140 * t142 - t146 * t236) * t128;
t119 = (t137 * t142 - t145 * t236) * t128;
t1 = [-t138 * t142 * t128, t120, 0, t119, 0, 0; (t241 * t122 - (-t126 + (-t142 * t238 + t126) * t128) * t235) * t121 ((-t157 * t215 + t158 * t202 + t169 * t212) * t122 - (t216 * t120 - t126 * t140 + t127 * t146) * t239) * t121, 0 (t139 * t122 - (t216 * t119 + t126 * t137 + t127 * t145) * t239) * t121, 0, 0; ((t137 * t201 - t147 * t205) * t130 - (t137 * t205 + t147 * t201) * t237) * t125 ((t141 * t201 - t149 * t205) * t130 - (t141 * t205 + t149 * t201) * t237) * t125, 0 (-t201 * t130 + t205 * t237) * t138 * t125, t225 * t125, 0;];
Ja_rot  = t1;
