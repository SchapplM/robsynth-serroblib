% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR10_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR10_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:13
% EndTime: 2019-02-26 22:21:14
% DurationCPUTime: 0.93s
% Computational Cost: add. (1265->118), mult. (4276->253), div. (493->12), fcn. (5073->11), ass. (0->110)
t180 = sin(qJ(2));
t172 = t180 ^ 2;
t184 = cos(qJ(2));
t175 = 0.1e1 / t184 ^ 2;
t237 = t172 * t175;
t181 = sin(qJ(1));
t173 = t181 ^ 2;
t167 = t173 * t237 + 0.1e1;
t174 = 0.1e1 / t184;
t234 = t174 * t180;
t256 = t180 * t237;
t193 = qJD(2) * (t174 * t256 + t234);
t185 = cos(qJ(1));
t226 = qJD(1) * t185;
t235 = t172 * t181;
t202 = t226 * t235;
t240 = (t173 * t193 + t175 * t202) / t167 ^ 2;
t257 = -0.2e1 * t240;
t209 = 0.1e1 + t237;
t255 = t181 * t209;
t183 = cos(qJ(3));
t254 = (-qJD(3) + qJD(5)) * t183;
t179 = sin(qJ(3));
t222 = qJD(3) * t179;
t253 = -qJD(5) * t179 + t222;
t229 = t181 * t184;
t194 = t179 * t229 + t183 * t185;
t223 = qJD(2) * t185;
t210 = t180 * t223;
t228 = t184 * t185;
t212 = t183 * t228;
t132 = t194 * qJD(1) - qJD(3) * t212 + t179 * t210 - t181 * t222;
t204 = -qJD(1) * t184 + qJD(3);
t205 = qJD(3) * t184 - qJD(1);
t233 = t179 * t185;
t133 = -t205 * t233 + (t204 * t181 - t210) * t183;
t178 = sin(qJ(5));
t182 = cos(qJ(5));
t230 = t181 * t183;
t160 = t179 * t228 - t230;
t161 = t181 * t179 + t212;
t198 = t160 * t182 - t161 * t178;
t125 = t198 * qJD(5) - t132 * t178 + t133 * t182;
t148 = t160 * t178 + t161 * t182;
t140 = 0.1e1 / t148;
t196 = t178 * t179 + t182 * t183;
t197 = t178 * t183 - t179 * t182;
t141 = 0.1e1 / t148 ^ 2;
t243 = t141 * t198;
t252 = t197 * t140 + t196 * t243;
t231 = t181 * t180;
t166 = atan2(t231, t184);
t163 = cos(t166);
t162 = sin(t166);
t214 = t162 * t231;
t152 = t163 * t184 + t214;
t149 = 0.1e1 / t152;
t150 = 0.1e1 / t152 ^ 2;
t251 = 0.2e1 * t180;
t164 = 0.1e1 / t167;
t250 = t164 - 0.1e1;
t124 = t148 * qJD(5) + t132 * t182 + t133 * t178;
t139 = t198 ^ 2;
t130 = t139 * t141 + 0.1e1;
t142 = t140 * t141;
t246 = t125 * t142;
t249 = (-t124 * t243 - t139 * t246) / t130 ^ 2;
t177 = t185 ^ 2;
t236 = t172 * t177;
t138 = t150 * t236 + 0.1e1;
t224 = qJD(2) * t184;
t211 = t180 * t226;
t225 = qJD(2) * t181;
t131 = ((t181 * t224 + t211) * t174 + t225 * t237) * t164;
t238 = t163 * t180;
t122 = (t131 * t181 - qJD(2)) * t238 + (t211 + (-t131 + t225) * t184) * t162;
t247 = t122 * t149 * t150;
t248 = (-t236 * t247 + (t177 * t180 * t224 - t202) * t150) / t138 ^ 2;
t245 = t131 * t162;
t244 = t131 * t180;
t232 = t180 * t185;
t156 = t196 * t232;
t242 = t141 * t156;
t241 = t150 * t180;
t154 = t164 * t255;
t239 = t154 * t181;
t227 = qJD(1) * t181;
t218 = 0.2e1 * t249;
t217 = -0.2e1 * t247;
t216 = -0.2e1 * t142 * t198;
t215 = t150 * t232;
t213 = t164 * t172 * t174;
t208 = -0.2e1 * t180 * t248;
t207 = t125 * t216;
t206 = t174 * t257;
t203 = t181 * t213;
t201 = t209 * t185;
t159 = -t183 * t229 + t233;
t199 = -t159 * t178 - t182 * t194;
t144 = t159 * t182 - t178 * t194;
t195 = t204 * t185;
t155 = t197 * t232;
t136 = 0.1e1 / t138;
t135 = t183 * t195 + (qJD(2) * t180 * t183 + t205 * t179) * t181;
t134 = -t205 * t230 + (t180 * t225 + t195) * t179;
t128 = 0.1e1 / t130;
t127 = (-t250 * t180 * t162 + t163 * t203) * t185;
t126 = t162 * t229 - t238 + (-t162 * t184 + t163 * t231) * t154;
t123 = t255 * t257 + (qJD(1) * t201 + 0.2e1 * t181 * t193) * t164;
t1 = [t206 * t232 + (qJD(2) * t201 - t227 * t234) * t164, t123, 0, 0, 0, 0; (t149 * t208 + (t149 * t224 + (-qJD(1) * t127 - t122) * t241) * t136) * t181 + (t150 * t208 * t127 + (((-t131 * t203 - t250 * t224 + t240 * t251) * t162 + (t206 * t235 + t244 + (-t244 + (t251 + t256) * t225) * t164) * t163) * t215 + (t150 * t224 + t180 * t217) * t127 + (t149 + ((-t173 + t177) * t163 * t213 + t250 * t214) * t150) * t180 * qJD(1)) * t136) * t185, 0.2e1 * (-t126 * t241 + t149 * t184) * t185 * t248 + ((t149 * t227 + (qJD(2) * t126 + t122) * t185 * t150) * t184 + (t149 * t223 + (t123 * t163 * t181 - t162 * t225 - t239 * t245 + t245 + (qJD(2) * t162 + t163 * t226) * t154) * t215 + (-t150 * t227 + t185 * t217) * t126 + ((-t123 + t226) * t162 + ((-0.1e1 + t239) * qJD(2) + (-t154 + t181) * t131) * t163) * t150 * t228) * t180) * t136, 0, 0, 0, 0; (t140 * t199 - t144 * t243) * t218 + ((t144 * qJD(5) - t134 * t182 + t135 * t178) * t140 + t144 * t207 + (t199 * t125 + (t199 * qJD(5) + t134 * t178 + t135 * t182) * t198 - t144 * t124) * t141) * t128 (t140 * t155 + t198 * t242) * t218 + (t124 * t242 + (t141 * t155 - t156 * t216) * t125 - t252 * t184 * t223 + (t252 * t227 + ((-t254 * t140 + t253 * t243) * t182 + (t253 * t140 + t254 * t243) * t178) * t185) * t180) * t128 (t140 * t148 + t198 * t243) * t218 + (-t125 * t140 - t198 * t207 + (0.2e1 * t198 * t124 + t148 * t125) * t141) * t128, 0, -0.2e1 * t249 - 0.2e1 * (t124 * t128 * t141 - (-t128 * t246 - t141 * t249) * t198) * t198, 0;];
JaD_rot  = t1;
