% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR12_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobia_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:08
% EndTime: 2019-02-26 21:21:10
% DurationCPUTime: 1.22s
% Computational Cost: add. (5345->78), mult. (15419->180), div. (115->9), fcn. (20897->21), ass. (0->95)
t187 = cos(pkin(6));
t180 = sin(pkin(14));
t234 = sin(qJ(1));
t214 = t234 * t180;
t184 = cos(pkin(14));
t196 = cos(qJ(1));
t217 = t196 * t184;
t176 = -t187 * t217 + t214;
t213 = t234 * t184;
t218 = t196 * t180;
t177 = t187 * t218 + t213;
t191 = sin(qJ(3));
t195 = cos(qJ(3));
t182 = sin(pkin(7));
t183 = sin(pkin(6));
t221 = t183 * t196;
t216 = t182 * t221;
t186 = cos(pkin(7));
t219 = t186 * t191;
t169 = t176 * t219 - t177 * t195 + t191 * t216;
t190 = sin(qJ(4));
t194 = cos(qJ(4));
t168 = t195 * (t176 * t186 + t216) + t177 * t191;
t174 = -t176 * t182 + t186 * t221;
t181 = sin(pkin(8));
t185 = cos(pkin(8));
t210 = t168 * t185 + t174 * t181;
t154 = t169 * t194 + t190 * t210;
t189 = sin(qJ(5));
t193 = cos(qJ(5));
t201 = t168 * t181 - t174 * t185;
t240 = t154 * t189 + t193 * t201;
t142 = t154 * t193 - t189 * t201;
t237 = t169 * t190 - t194 * t210;
t222 = t182 * t187;
t173 = t191 * t222 + (t180 * t195 + t184 * t219) * t183;
t172 = t195 * t222 + (t184 * t186 * t195 - t180 * t191) * t183;
t175 = -t182 * t183 * t184 + t186 * t187;
t209 = t172 * t185 + t175 * t181;
t161 = t173 * t194 + t190 * t209;
t166 = -t172 * t181 + t175 * t185;
t149 = t161 * t189 - t166 * t193;
t138 = atan2(t240, t149);
t133 = sin(t138);
t134 = cos(t138);
t129 = t133 * t240 + t134 * t149;
t128 = 0.1e1 / t129 ^ 2;
t207 = -t187 * t213 - t218;
t215 = t183 * t234;
t203 = t182 * t215 + t186 * t207;
t206 = -t187 * t214 + t217;
t199 = t191 * t206 - t195 * t203;
t204 = -t182 * t207 + t186 * t215;
t202 = t204 * t181;
t170 = t191 * t203 + t195 * t206;
t227 = t170 * t194;
t156 = t227 + (-t185 * t199 + t202) * t190;
t197 = t181 * t199 + t185 * t204;
t143 = t156 * t189 - t193 * t197;
t233 = t128 * t143;
t144 = t156 * t193 + t189 * t197;
t198 = t199 * t194;
t155 = t170 * t190 + t185 * t198 - t194 * t202;
t188 = sin(qJ(6));
t192 = cos(qJ(6));
t136 = t144 * t192 + t155 * t188;
t132 = 0.1e1 / t136 ^ 2;
t135 = t144 * t188 - t155 * t192;
t232 = t132 * t135;
t231 = t134 * t240;
t148 = 0.1e1 / t149 ^ 2;
t230 = t240 * t148;
t229 = t143 ^ 2 * t128;
t228 = t155 * t193;
t223 = t181 * t193;
t220 = t185 * t190;
t212 = t132 * t135 ^ 2 + 0.1e1;
t211 = -t133 * t149 + t231;
t160 = -t173 * t190 + t194 * t209;
t159 = (t172 * t194 - t173 * t220) * t189 - t173 * t223;
t158 = -t170 * t220 - t198;
t157 = t185 * t227 - t190 * t199;
t150 = t161 * t193 + t166 * t189;
t147 = 0.1e1 / t149;
t146 = t170 * t181 * t189 + t158 * t193;
t145 = (-t168 * t194 + t169 * t220) * t189 + t169 * t223;
t137 = 0.1e1 / (t148 * t240 ^ 2 + 0.1e1);
t131 = 0.1e1 / t136;
t130 = 0.1e1 / t212;
t127 = 0.1e1 / t129;
t126 = 0.1e1 / (0.1e1 + t229);
t125 = (-t147 * t237 - t160 * t230) * t189 * t137;
t124 = (-t145 * t147 - t159 * t230) * t137;
t123 = (t142 * t147 - t150 * t230) * t137;
t1 = [-t143 * t147 * t137, 0, t124, t125, t123, 0; (t240 * t127 - (-t133 + (-t147 * t231 + t133) * t137) * t229) * t126, 0 ((t158 * t189 - t170 * t223) * t127 - (t124 * t211 - t133 * t145 + t134 * t159) * t233) * t126 (-t155 * t189 * t127 - ((-t133 * t237 + t134 * t160) * t189 + t211 * t125) * t233) * t126 (t144 * t127 - (t123 * t211 + t133 * t142 + t134 * t150) * t233) * t126, 0; ((t142 * t188 - t192 * t237) * t131 - (t142 * t192 + t188 * t237) * t232) * t130, 0 ((t146 * t188 - t157 * t192) * t131 - (t146 * t192 + t157 * t188) * t232) * t130 ((-t156 * t192 - t188 * t228) * t131 - (t156 * t188 - t192 * t228) * t232) * t130 (-t131 * t188 + t192 * t232) * t143 * t130, t212 * t130;];
Ja_rot  = t1;
