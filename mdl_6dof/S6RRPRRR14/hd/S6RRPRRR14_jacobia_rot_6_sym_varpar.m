% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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
% Datum: 2019-01-03 10:25
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR14_jacobia_rot_6_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_rot_6_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_rot_6_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-03 10:25:34
% EndTime: 2019-01-03 10:25:35
% DurationCPUTime: 1.40s
% Computational Cost: add. (5348->88), mult. (15437->209), div. (115->9), fcn. (20918->21), ass. (0->100)
t203 = sin(qJ(2));
t207 = cos(qJ(2));
t208 = cos(qJ(1));
t247 = cos(pkin(6));
t226 = t208 * t247;
t248 = sin(qJ(1));
t189 = t203 * t226 + t248 * t207;
t193 = sin(pkin(14));
t197 = cos(pkin(14));
t188 = t248 * t203 - t207 * t226;
t195 = sin(pkin(7));
t199 = cos(pkin(7));
t196 = sin(pkin(6));
t233 = t196 * t208;
t219 = t188 * t199 + t195 * t233;
t177 = -t189 * t197 + t219 * t193;
t202 = sin(qJ(4));
t206 = cos(qJ(4));
t176 = t189 * t193 + t219 * t197;
t185 = -t188 * t195 + t199 * t233;
t194 = sin(pkin(8));
t198 = cos(pkin(8));
t222 = t176 * t198 + t185 * t194;
t161 = t177 * t206 + t222 * t202;
t201 = sin(qJ(5));
t205 = cos(qJ(5));
t212 = t176 * t194 - t185 * t198;
t255 = t161 * t201 + t212 * t205;
t149 = t161 * t205 - t212 * t201;
t252 = t177 * t202 - t222 * t206;
t224 = t247 * t248;
t190 = -t203 * t224 + t208 * t207;
t218 = t208 * t203 + t207 * t224;
t228 = t196 * t248;
t215 = t195 * t228 - t218 * t199;
t213 = t190 * t193 - t215 * t197;
t216 = t218 * t195 + t199 * t228;
t249 = -t216 * t194 + t213 * t198;
t225 = t247 * t195;
t231 = t199 * t207;
t184 = t196 * t203 * t197 + (t196 * t231 + t225) * t193;
t183 = t197 * t225 + (-t193 * t203 + t197 * t231) * t196;
t187 = -t196 * t207 * t195 + t247 * t199;
t221 = t183 * t198 + t187 * t194;
t168 = t184 * t206 + t221 * t202;
t174 = -t183 * t194 + t187 * t198;
t156 = t168 * t201 - t174 * t205;
t143 = atan2(t255, t156);
t138 = sin(t143);
t139 = cos(t143);
t136 = t138 * t255 + t139 * t156;
t135 = 0.1e1 / t136 ^ 2;
t178 = t190 * t197 + t215 * t193;
t163 = t178 * t206 - t249 * t202;
t209 = t213 * t194 + t216 * t198;
t150 = t163 * t201 - t209 * t205;
t246 = t135 * t150;
t245 = t139 * t255;
t151 = t163 * t205 + t209 * t201;
t162 = t178 * t202 + t249 * t206;
t200 = sin(qJ(6));
t204 = cos(qJ(6));
t145 = t151 * t204 + t162 * t200;
t142 = 0.1e1 / t145 ^ 2;
t144 = t151 * t200 - t162 * t204;
t244 = t142 * t144;
t155 = 0.1e1 / t156 ^ 2;
t243 = t255 * t155;
t242 = t150 ^ 2 * t135;
t241 = t162 * t205;
t236 = t193 * t199;
t235 = t195 * t194;
t234 = t195 * t198;
t232 = t197 * t199;
t230 = t203 * t195;
t229 = t203 * t199;
t227 = t144 ^ 2 * t142 + 0.1e1;
t223 = -t138 * t156 + t245;
t180 = -t190 * t232 + t218 * t193;
t220 = t180 * t198 + t190 * t235;
t181 = -t190 * t236 - t218 * t197;
t179 = t188 * t193 - t189 * t232;
t172 = -t180 * t194 + t190 * t234;
t167 = -t184 * t202 + t221 * t206;
t166 = ((t198 * t202 * t201 + t194 * t205) * (-t193 * t207 - t197 * t229) + ((-t193 * t229 + t197 * t207) * t206 + t194 * t202 * t230) * t201 - t198 * t205 * t230) * t196;
t165 = t181 * t206 + t220 * t202;
t164 = t181 * t202 - t220 * t206;
t157 = t168 * t205 + t174 * t201;
t154 = 0.1e1 / t156;
t153 = t165 * t205 + t172 * t201;
t152 = ((-t188 * t197 - t189 * t236) * t206 + (t179 * t198 + t189 * t235) * t202) * t201 - (-t179 * t194 + t189 * t234) * t205;
t141 = 0.1e1 / t145;
t140 = 0.1e1 / (t155 * t255 ^ 2 + 0.1e1);
t137 = 0.1e1 / t227;
t134 = 0.1e1 / t136;
t133 = 0.1e1 / (0.1e1 + t242);
t132 = (-t154 * t252 - t167 * t243) * t201 * t140;
t131 = (-t152 * t154 - t166 * t243) * t140;
t130 = (t149 * t154 - t157 * t243) * t140;
t1 = [-t150 * t154 * t140, t131, 0, t132, t130, 0; (t255 * t134 - (-t138 + (-t154 * t245 + t138) * t140) * t242) * t133 ((t165 * t201 - t172 * t205) * t134 - (t223 * t131 - t138 * t152 + t139 * t166) * t246) * t133, 0 (-t162 * t201 * t134 - ((-t138 * t252 + t139 * t167) * t201 + t223 * t132) * t246) * t133 (t151 * t134 - (t223 * t130 + t138 * t149 + t139 * t157) * t246) * t133, 0; ((t149 * t200 - t204 * t252) * t141 - (t149 * t204 + t200 * t252) * t244) * t137 ((t153 * t200 - t164 * t204) * t141 - (t153 * t204 + t164 * t200) * t244) * t137, 0 ((-t163 * t204 - t200 * t241) * t141 - (t163 * t200 - t204 * t241) * t244) * t137 (-t200 * t141 + t204 * t244) * t150 * t137, t227 * t137;];
Ja_rot  = t1;
