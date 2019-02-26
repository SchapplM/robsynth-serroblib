% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 7 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% Ja_rot [3x7]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S7RRRRRRR1_jacobia_rot_7_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobia_rot_7_floatb_twist_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobia_rot_7_floatb_twist_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_7_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 21:20:54
% EndTime: 2018-11-26 21:20:55
% DurationCPUTime: 1.04s
% Computational Cost: add. (2869->90), mult. (7788->218), div. (175->9), fcn. (10858->17), ass. (0->100)
t196 = sin(qJ(1));
t200 = cos(qJ(4));
t201 = cos(qJ(3));
t194 = sin(qJ(3));
t237 = cos(qJ(1));
t214 = t237 * t194;
t202 = cos(qJ(2));
t221 = t196 * t202;
t207 = t201 * t221 + t214;
t193 = sin(qJ(4));
t195 = sin(qJ(2));
t224 = t195 * t193;
t175 = t196 * t224 + t200 * t207;
t213 = t237 * t201;
t184 = t194 * t221 - t213;
t192 = sin(qJ(5));
t199 = cos(qJ(5));
t159 = t175 * t199 - t184 * t192;
t191 = sin(qJ(6));
t198 = cos(qJ(6));
t223 = t195 * t200;
t203 = t193 * t207 - t196 * t223;
t142 = t159 * t191 - t198 * t203;
t143 = t159 * t198 + t191 * t203;
t238 = -t175 * t192 - t184 * t199;
t219 = t202 * t193;
t222 = t195 * t201;
t183 = t200 * t222 - t219;
t225 = t194 * t195;
t174 = t183 * t199 - t192 * t225;
t218 = t202 * t200;
t182 = t193 * t222 + t218;
t156 = -t174 * t191 + t182 * t198;
t141 = atan2(t142, t156);
t138 = sin(t141);
t139 = cos(t141);
t132 = t138 * t142 + t139 * t156;
t131 = 0.1e1 / t132 ^ 2;
t208 = -t196 * t194 + t202 * t213;
t216 = t195 * t237;
t210 = t193 * t216;
t178 = t200 * t208 + t210;
t185 = -t196 * t201 - t202 * t214;
t228 = t185 * t192;
t164 = t178 * t199 + t228;
t187 = t200 * t216;
t206 = t193 * t208 - t187;
t145 = t164 * t191 - t198 * t206;
t236 = t131 * t145;
t147 = t164 * t198 + t191 * t206;
t163 = t178 * t192 - t185 * t199;
t190 = sin(qJ(7));
t197 = cos(qJ(7));
t137 = t147 * t197 - t163 * t190;
t135 = 0.1e1 / t137 ^ 2;
t136 = t147 * t190 + t163 * t197;
t235 = t135 * t136;
t234 = t139 * t142;
t155 = 0.1e1 / t156 ^ 2;
t233 = t142 * t155;
t232 = t145 ^ 2 * t131;
t231 = t163 * t198;
t227 = t191 * t199;
t226 = t193 * t198;
t220 = t200 * t199;
t215 = t202 * t237;
t212 = t136 ^ 2 * t135 + 0.1e1;
t211 = t195 * t214;
t209 = -t138 * t156 + t234;
t205 = t192 * t206;
t204 = t199 * t206;
t180 = -t187 * t201 + t193 * t215;
t179 = -t200 * t215 - t201 * t210;
t173 = -t183 * t192 - t199 * t225;
t170 = t180 * t199 + t192 * t211;
t169 = t180 * t192 - t199 * t211;
t168 = (-(-t192 * t201 - t194 * t220) * t191 - t194 * t226) * t195;
t167 = t185 * t220 - t192 * t208;
t166 = t199 * t208 + t200 * t228;
t165 = t182 * t227 + t183 * t198;
t162 = -((t201 * t218 + t224) * t199 - t202 * t194 * t192) * t191 + (t201 * t219 - t223) * t198;
t157 = -t174 * t198 - t182 * t191;
t154 = 0.1e1 / t156;
t153 = t170 * t198 + t179 * t191;
t152 = t156 * t196;
t151 = t185 * t193 * t191 + t167 * t198;
t150 = (-t184 * t220 - t192 * t207) * t191 + t184 * t226;
t149 = t178 * t191 - t198 * t204;
t148 = -t175 * t198 - t203 * t227;
t140 = 0.1e1 / (t142 ^ 2 * t155 + 0.1e1);
t134 = 0.1e1 / t137;
t133 = 0.1e1 / t212;
t130 = 0.1e1 / t132;
t129 = 0.1e1 / (0.1e1 + t232);
t128 = (t154 * t238 + t173 * t233) * t191 * t140;
t127 = (t150 * t154 - t168 * t233) * t140;
t126 = (t148 * t154 - t165 * t233) * t140;
t125 = (t152 * t154 - t162 * t233) * t140;
t124 = (t143 * t154 - t157 * t233) * t140;
t1 = [t145 * t154 * t140, t125, t127, t126, t128, t124, 0; (t142 * t130 + (t138 + (t154 * t234 - t138) * t140) * t232) * t129 ((-t170 * t191 + t179 * t198) * t130 + (t125 * t209 + t138 * t152 + t139 * t162) * t236) * t129 ((-t167 * t191 + t185 * t226) * t130 + (t127 * t209 + t138 * t150 + t139 * t168) * t236) * t129 ((t178 * t198 + t191 * t204) * t130 + (t126 * t209 + t138 * t148 + t139 * t165) * t236) * t129 (t163 * t191 * t130 + ((t138 * t238 - t139 * t173) * t191 + t209 * t128) * t236) * t129 (-t147 * t130 + (t124 * t209 + t138 * t143 + t139 * t157) * t236) * t129, 0; ((-t143 * t190 + t197 * t238) * t134 - (-t143 * t197 - t190 * t238) * t235) * t133 ((t153 * t190 + t169 * t197) * t134 - (t153 * t197 - t169 * t190) * t235) * t133 ((t151 * t190 + t166 * t197) * t134 - (t151 * t197 - t166 * t190) * t235) * t133 ((t149 * t190 - t197 * t205) * t134 - (t149 * t197 + t190 * t205) * t235) * t133 ((t164 * t197 - t190 * t231) * t134 - (-t164 * t190 - t197 * t231) * t235) * t133 (-t190 * t134 + t197 * t235) * t145 * t133, t212 * t133;];
Ja_rot  = t1;
