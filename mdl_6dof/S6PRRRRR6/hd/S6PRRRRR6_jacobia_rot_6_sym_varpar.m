% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRR6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobia_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:52
% EndTime: 2019-02-26 20:21:52
% DurationCPUTime: 0.87s
% Computational Cost: add. (5587->99), mult. (16161->229), div. (125->9), fcn. (21882->21), ass. (0->109)
t187 = sin(pkin(14));
t191 = cos(pkin(14));
t199 = sin(qJ(2));
t194 = cos(pkin(6));
t204 = cos(qJ(2));
t222 = t194 * t204;
t181 = -t187 * t199 + t191 * t222;
t223 = t194 * t199;
t182 = t187 * t204 + t191 * t223;
t198 = sin(qJ(3));
t203 = cos(qJ(3));
t189 = sin(pkin(7));
t190 = sin(pkin(6));
t232 = t189 * t190;
t217 = t191 * t232;
t193 = cos(pkin(7));
t225 = t193 * t198;
t170 = t181 * t225 + t182 * t203 - t198 * t217;
t197 = sin(qJ(4));
t202 = cos(qJ(4));
t188 = sin(pkin(8));
t192 = cos(pkin(8));
t207 = -t182 * t198 + (t181 * t193 - t217) * t203;
t228 = t190 * t193;
t212 = -t181 * t189 - t191 * t228;
t205 = t188 * t212 + t192 * t207;
t154 = t170 * t202 + t197 * t205;
t196 = sin(qJ(5));
t201 = cos(qJ(5));
t206 = -t188 * t207 + t192 * t212;
t140 = t154 * t196 - t201 * t206;
t219 = t199 * t203;
t220 = t198 * t204;
t230 = t189 * t194;
t178 = t198 * t230 + (t193 * t220 + t219) * t190;
t218 = t203 * t204;
t221 = t198 * t199;
t177 = t203 * t230 + (t193 * t218 - t221) * t190;
t180 = t194 * t193 - t204 * t232;
t214 = t177 * t192 + t180 * t188;
t166 = t178 * t202 + t197 * t214;
t169 = -t177 * t188 + t180 * t192;
t151 = t166 * t196 - t169 * t201;
t139 = atan2(-t140, t151);
t136 = sin(t139);
t137 = cos(t139);
t130 = -t136 * t140 + t137 * t151;
t129 = 0.1e1 / t130 ^ 2;
t184 = -t187 * t223 + t191 * t204;
t183 = -t187 * t222 - t191 * t199;
t210 = t183 * t193 + t187 * t232;
t171 = -t184 * t198 + t203 * t210;
t172 = t184 * t203 + t198 * t210;
t211 = -t183 * t189 + t187 * t228;
t209 = t211 * t188;
t156 = t172 * t202 + (t171 * t192 + t209) * t197;
t208 = -t171 * t188 + t192 * t211;
t143 = t156 * t196 - t201 * t208;
t238 = t129 * t143;
t144 = t156 * t201 + t196 * t208;
t226 = t192 * t202;
t155 = -t171 * t226 + t172 * t197 - t202 * t209;
t195 = sin(qJ(6));
t200 = cos(qJ(6));
t135 = t144 * t200 + t155 * t195;
t133 = 0.1e1 / t135 ^ 2;
t134 = t144 * t195 - t155 * t200;
t237 = t133 * t134;
t150 = 0.1e1 / t151 ^ 2;
t236 = t140 * t150;
t235 = t155 * t201;
t234 = t188 * t201;
t233 = t189 * t188;
t231 = t189 * t192;
t229 = t189 * t199;
t227 = t192 * t197;
t224 = t193 * t203;
t216 = t134 ^ 2 * t133 + 0.1e1;
t215 = -t136 * t151 - t137 * t140;
t174 = -t183 * t198 - t184 * t224;
t213 = t174 * t192 + t184 * t233;
t175 = t183 * t203 - t184 * t225;
t173 = -t181 * t198 - t182 * t224;
t167 = -t174 * t188 + t184 * t231;
t165 = -t178 * t197 + t202 * t214;
t162 = (t177 * t202 - t178 * t227) * t196 - t178 * t234;
t161 = ((t196 * t227 + t234) * (-t193 * t219 - t220) + ((-t193 * t221 + t218) * t202 + t188 * t197 * t229) * t196 - t192 * t201 * t229) * t190;
t160 = t171 * t202 - t172 * t227;
t159 = t171 * t197 + t172 * t226;
t158 = t175 * t202 + t197 * t213;
t157 = t175 * t197 - t202 * t213;
t153 = -t170 * t197 + t202 * t205;
t152 = t166 * t201 + t169 * t196;
t149 = 0.1e1 / t151;
t148 = t172 * t188 * t196 + t160 * t201;
t147 = (-t170 * t227 + t202 * t207) * t196 - t170 * t234;
t146 = t158 * t201 + t167 * t196;
t145 = ((t181 * t203 - t182 * t225) * t202 + (t173 * t192 + t182 * t233) * t197) * t196 - (-t173 * t188 + t182 * t231) * t201;
t142 = t154 * t201 + t196 * t206;
t138 = 0.1e1 / (t140 ^ 2 * t150 + 0.1e1);
t132 = 0.1e1 / t135;
t131 = 0.1e1 / t216;
t128 = 0.1e1 / t130;
t127 = 0.1e1 / (t143 ^ 2 * t129 + 0.1e1);
t126 = (-t149 * t153 + t165 * t236) * t196 * t138;
t125 = (-t147 * t149 + t162 * t236) * t138;
t124 = (-t145 * t149 + t161 * t236) * t138;
t123 = (-t142 * t149 + t152 * t236) * t138;
t1 = [0, t124, t125, t126, t123, 0; 0 ((t158 * t196 - t167 * t201) * t128 - (t124 * t215 - t136 * t145 + t137 * t161) * t238) * t127 ((t160 * t196 - t172 * t234) * t128 - (t125 * t215 - t136 * t147 + t137 * t162) * t238) * t127 (-t155 * t196 * t128 - ((-t136 * t153 + t137 * t165) * t196 + t215 * t126) * t238) * t127 (t144 * t128 - (t123 * t215 - t136 * t142 + t137 * t152) * t238) * t127, 0; 0 ((t146 * t195 - t157 * t200) * t132 - (t146 * t200 + t157 * t195) * t237) * t131 ((t148 * t195 - t159 * t200) * t132 - (t148 * t200 + t159 * t195) * t237) * t131 ((-t156 * t200 - t195 * t235) * t132 - (t156 * t195 - t200 * t235) * t237) * t131 (-t195 * t132 + t200 * t237) * t143 * t131, t216 * t131;];
Ja_rot  = t1;
