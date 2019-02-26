% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:02
% EndTime: 2019-02-26 21:02:03
% DurationCPUTime: 1.01s
% Computational Cost: add. (2901->118), mult. (4276->251), div. (493->12), fcn. (5073->11), ass. (0->114)
t188 = sin(qJ(3));
t182 = t188 ^ 2;
t191 = cos(qJ(3));
t184 = 0.1e1 / t191 ^ 2;
t239 = t182 * t184;
t180 = qJ(1) + pkin(10);
t178 = sin(t180);
t176 = t178 ^ 2;
t172 = t176 * t239 + 0.1e1;
t183 = 0.1e1 / t191;
t265 = t188 * t239;
t199 = qJD(3) * (t188 + t265) * t183;
t179 = cos(t180);
t235 = qJD(1) * t179;
t243 = t178 * t182;
t208 = t235 * t243;
t250 = (t176 * t199 + t184 * t208) / t172 ^ 2;
t267 = -0.2e1 * t250;
t215 = 0.1e1 + t239;
t266 = t178 * t215;
t190 = cos(qJ(4));
t264 = (-qJD(4) + qJD(6)) * t190;
t187 = sin(qJ(4));
t230 = qJD(4) * t187;
t263 = -qJD(6) * t187 + t230;
t211 = qJD(4) * t191 - qJD(1);
t232 = qJD(3) * t188;
t262 = t187 * t211 + t190 * t232;
t238 = t187 * t191;
t200 = t178 * t238 + t179 * t190;
t217 = t187 * t232;
t237 = t190 * t191;
t220 = t179 * t237;
t137 = qJD(1) * t200 - qJD(4) * t220 - t178 * t230 + t179 * t217;
t210 = -qJD(1) * t191 + qJD(4);
t201 = t210 * t190;
t138 = t178 * t201 - t179 * t262;
t186 = sin(qJ(6));
t189 = cos(qJ(6));
t165 = -t178 * t190 + t179 * t238;
t166 = t178 * t187 + t220;
t204 = t165 * t189 - t166 * t186;
t130 = qJD(6) * t204 - t137 * t186 + t138 * t189;
t153 = t165 * t186 + t166 * t189;
t145 = 0.1e1 / t153;
t202 = t186 * t187 + t189 * t190;
t203 = t186 * t190 - t187 * t189;
t146 = 0.1e1 / t153 ^ 2;
t252 = t146 * t204;
t261 = t145 * t203 + t202 * t252;
t242 = t178 * t188;
t171 = atan2(t242, t191);
t168 = cos(t171);
t167 = sin(t171);
t222 = t167 * t242;
t159 = t168 * t191 + t222;
t156 = 0.1e1 / t159;
t157 = 0.1e1 / t159 ^ 2;
t260 = 0.2e1 * t188;
t169 = 0.1e1 / t172;
t259 = t169 - 0.1e1;
t129 = qJD(6) * t153 + t137 * t189 + t138 * t186;
t144 = t204 ^ 2;
t135 = t144 * t146 + 0.1e1;
t147 = t145 * t146;
t255 = t130 * t147;
t258 = (-t129 * t252 - t144 * t255) / t135 ^ 2;
t177 = t179 ^ 2;
t244 = t177 * t182;
t143 = t157 * t244 + 0.1e1;
t231 = qJD(3) * t191;
t234 = qJD(1) * t188;
t218 = t179 * t234;
t233 = qJD(3) * t178;
t136 = ((t178 * t231 + t218) * t183 + t233 * t239) * t169;
t245 = t168 * t188;
t127 = (t136 * t178 - qJD(3)) * t245 + (t218 + (-t136 + t233) * t191) * t167;
t256 = t127 * t156 * t157;
t257 = (-t244 * t256 + (t177 * t188 * t231 - t208) * t157) / t143 ^ 2;
t254 = t136 * t167;
t253 = t136 * t188;
t240 = t179 * t188;
t161 = t202 * t240;
t251 = t146 * t161;
t155 = t169 * t266;
t249 = t155 * t178;
t248 = t157 * t179;
t247 = t157 * t188;
t246 = t167 * t191;
t241 = t179 * t187;
t236 = qJD(1) * t178;
t226 = 0.2e1 * t258;
t225 = -0.2e1 * t256;
t224 = -0.2e1 * t147 * t204;
t223 = t157 * t240;
t221 = t169 * t182 * t183;
t219 = t178 * t234;
t214 = -0.2e1 * t188 * t257;
t213 = t130 * t224;
t212 = t183 * t267;
t209 = t178 * t221;
t207 = t215 * t179;
t164 = -t178 * t237 + t241;
t205 = -t164 * t186 - t189 * t200;
t149 = t164 * t189 - t186 * t200;
t160 = t203 * t240;
t141 = 0.1e1 / t143;
t140 = t178 * t262 + t179 * t201;
t139 = t210 * t241 + (-t190 * t211 + t217) * t178;
t133 = 0.1e1 / t135;
t132 = (-t167 * t188 * t259 + t168 * t209) * t179;
t131 = t178 * t246 - t245 + (t168 * t242 - t246) * t155;
t128 = t266 * t267 + (qJD(1) * t207 + 0.2e1 * t178 * t199) * t169;
t1 = [t212 * t240 + (qJD(3) * t207 - t183 * t219) * t169, 0, t128, 0, 0, 0; (t156 * t214 + (t156 * t231 + (-qJD(1) * t132 - t127) * t247) * t141) * t178 + (t157 * t214 * t132 + (((-t136 * t209 - t231 * t259 + t250 * t260) * t167 + (t212 * t243 + t253 + (-t253 + (t260 + t265) * t233) * t169) * t168) * t223 + (t157 * t231 + t188 * t225) * t132 + (t156 + ((-t176 + t177) * t168 * t221 + t259 * t222) * t157) * t234) * t141) * t179, 0, 0.2e1 * (-t131 * t247 + t156 * t191) * t179 * t257 + ((t156 * t236 + (qJD(3) * t131 + t127) * t248) * t191 + (t179 * qJD(3) * t156 + (t128 * t168 * t178 - t167 * t233 - t249 * t254 + t254 + (qJD(3) * t167 + t168 * t235) * t155) * t223 + (-t157 * t236 + t179 * t225) * t131 + ((-t128 + t235) * t167 + ((-0.1e1 + t249) * qJD(3) + (-t155 + t178) * t136) * t168) * t191 * t248) * t188) * t141, 0, 0, 0; (t145 * t205 - t149 * t252) * t226 + ((qJD(6) * t149 - t139 * t189 + t140 * t186) * t145 + t149 * t213 + (t205 * t130 + (qJD(6) * t205 + t139 * t186 + t140 * t189) * t204 - t149 * t129) * t146) * t133, 0 (t145 * t160 + t204 * t251) * t226 + (t129 * t251 + (t146 * t160 - t161 * t224) * t130 + t261 * t219 + (-t261 * t231 + ((-t145 * t264 + t252 * t263) * t189 + (t145 * t263 + t252 * t264) * t186) * t188) * t179) * t133 (t145 * t153 + t204 * t252) * t226 + (-t130 * t145 - t204 * t213 + (0.2e1 * t129 * t204 + t130 * t153) * t146) * t133, 0, -0.2e1 * t258 - 0.2e1 * (t129 * t133 * t146 - (-t133 * t255 - t146 * t258) * t204) * t204;];
JaD_rot  = t1;
