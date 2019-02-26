% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:58:03
% EndTime: 2019-02-26 20:58:05
% DurationCPUTime: 1.09s
% Computational Cost: add. (6874->125), mult. (6168->272), div. (1114->15), fcn. (7752->9), ass. (0->116)
t176 = pkin(9) + qJ(3);
t174 = cos(t176);
t177 = qJ(4) + pkin(10);
t173 = sin(t177);
t251 = sin(qJ(1));
t211 = t251 * t173;
t175 = cos(t177);
t179 = cos(qJ(1));
t231 = t179 * t175;
t154 = t174 * t231 + t211;
t148 = 0.1e1 / t154 ^ 2;
t172 = sin(t176);
t165 = t172 ^ 2;
t178 = t179 ^ 2;
t239 = t165 * t178;
t216 = t148 * t239;
t144 = 0.1e1 + t216;
t203 = qJD(1) * t251;
t228 = qJD(3) * t179;
t207 = t172 * t228;
t189 = t174 * t203 + t207;
t202 = t251 * qJD(4);
t232 = t179 * t173;
t133 = (-qJD(4) * t174 + qJD(1)) * t232 + (t202 - t189) * t175;
t147 = 0.1e1 / t154;
t246 = t133 * t147 * t148;
t197 = t239 * t246;
t208 = qJD(3) * t172 * t178;
t254 = (-t197 + (-t165 * t179 * t203 + t174 * t208) * t148) / t144 ^ 2;
t233 = t172 * t179;
t150 = t174 * t211 + t231;
t194 = t173 * t202;
t225 = qJD(4) * t179;
t205 = t175 * t225;
t132 = t150 * qJD(1) + t173 * t207 - t174 * t205 - t194;
t210 = t251 * t175;
t153 = t174 * t232 - t210;
t166 = 0.1e1 / t172;
t169 = 0.1e1 / t173;
t170 = 0.1e1 / t173 ^ 2;
t226 = qJD(4) * t175;
t206 = t170 * t226;
t167 = 0.1e1 / t172 ^ 2;
t229 = qJD(3) * t174;
t209 = t167 * t229;
t238 = t166 * t169;
t253 = (t166 * t206 + t169 * t209) * t153 + t132 * t238;
t234 = t172 * t173;
t140 = atan2(-t150, t234);
t137 = cos(t140);
t136 = sin(t140);
t245 = t136 * t150;
t131 = t137 * t234 - t245;
t128 = 0.1e1 / t131;
t129 = 0.1e1 / t131 ^ 2;
t252 = 0.2e1 * t153;
t145 = t150 ^ 2;
t237 = t167 * t170;
t141 = t145 * t237 + 0.1e1;
t138 = 0.1e1 / t141;
t190 = t172 * t226 + t173 * t229;
t214 = t150 * t237;
t212 = t172 * t251;
t195 = qJD(3) * t212;
t196 = t175 * t203;
t230 = qJD(1) * t179;
t134 = t175 * t202 * t174 - t196 + (t230 * t174 - t195 - t225) * t173;
t217 = t134 * t238;
t120 = (t190 * t214 - t217) * t138;
t187 = -t120 * t150 + t190;
t116 = (-t120 * t234 - t134) * t136 + t187 * t137;
t130 = t128 * t129;
t250 = t116 * t130;
t168 = t166 / t165;
t171 = t169 * t170;
t249 = (t134 * t214 + (-t167 * t171 * t226 - t168 * t170 * t229) * t145) / t141 ^ 2;
t248 = t129 * t153;
t247 = t132 * t129;
t244 = t136 * t153;
t243 = t136 * t172;
t242 = t137 * t150;
t241 = t137 * t153;
t240 = t137 * t174;
t236 = t167 * t174;
t235 = t170 * t175;
t227 = qJD(4) * t173;
t146 = t153 ^ 2;
t126 = t129 * t146 + 0.1e1;
t224 = 0.2e1 * (-t146 * t250 - t153 * t247) / t126 ^ 2;
t223 = -0.2e1 * t249;
t222 = 0.2e1 * t254;
t221 = t130 * t252;
t220 = t166 * t249;
t219 = t129 * t244;
t215 = t150 * t238;
t213 = t169 * t236;
t192 = t150 * t213 + t251;
t127 = t192 * t138;
t204 = t251 - t127;
t201 = t128 * t224;
t200 = t129 * t224;
t199 = t233 * t252;
t198 = t169 * t220;
t152 = t174 * t210 - t232;
t193 = t150 * t235 - t152 * t169;
t191 = t148 * t152 * t179 - t251 * t147;
t142 = 0.1e1 / t144;
t135 = t154 * qJD(1) - t174 * t194 - t175 * t195 - t205;
t124 = 0.1e1 / t126;
t123 = t193 * t166 * t138;
t119 = (-t136 + (t137 * t215 + t136) * t138) * t153;
t118 = -t127 * t242 + (t204 * t243 + t240) * t173;
t117 = t137 * t172 * t175 - t136 * t152 + (-t136 * t234 - t242) * t123;
t115 = t192 * t223 + (t134 * t213 + t230 + (-t206 * t236 + (-0.2e1 * t168 * t174 ^ 2 - t166) * t169 * qJD(3)) * t150) * t138;
t113 = -0.2e1 * t193 * t220 + (-t193 * t209 + (t134 * t235 - t135 * t169 + (t152 * t235 + (-0.2e1 * t171 * t175 ^ 2 - t169) * t150) * qJD(4)) * t166) * t138;
t1 = [t253 * t138 + t198 * t252, 0, t115, t113, 0, 0; t150 * t201 + (-t134 * t128 + (t116 * t150 + t119 * t132) * t129) * t124 + (t119 * t200 + (0.2e1 * t119 * t250 + (t132 * t138 - t132 - (-t120 * t138 * t215 + t223) * t153) * t129 * t136 + (-(-0.2e1 * t150 * t198 - t120) * t248 + (-(t120 + t217) * t153 + t253 * t150) * t129 * t138) * t137) * t124) * t153, 0, t118 * t153 * t200 + (-(-t115 * t242 + (t120 * t245 - t134 * t137) * t127) * t248 + (t116 * t221 + t247) * t118 + (-t128 * t233 - (-t127 * t243 + t136 * t212 + t240) * t248) * t226) * t124 + (t201 * t233 + ((-t128 * t228 - (t204 * qJD(3) - t120) * t219) * t174 + (t128 * t203 + (t179 * t116 - (-t115 + t230) * t244 - (t204 * t120 - qJD(3)) * t241) * t129) * t172) * t124) * t173 (t117 * t248 - t128 * t154) * t224 + (t117 * t247 + t133 * t128 + (t117 * t221 - t154 * t129) * t116 - (t175 * t229 - t172 * t227 - t113 * t150 - t123 * t134 + (-t123 * t234 - t152) * t120) * t129 * t241 - (-t135 + (-t113 * t173 - t120 * t175) * t172 - t187 * t123) * t219) * t124, 0, 0; t191 * t172 * t222 + (-t191 * t229 + ((qJD(1) * t147 + 0.2e1 * t152 * t246) * t179 + (-t251 * t133 - t135 * t179 + t152 * t203) * t148) * t172) * t142, 0 (t147 * t174 * t179 + t175 * t216) * t222 + (0.2e1 * t175 * t197 + t189 * t147 + ((t133 * t179 - 0.2e1 * t175 * t208) * t174 + (t178 * t227 + 0.2e1 * t179 * t196) * t165) * t148) * t142, t148 * t199 * t254 + (t199 * t246 + (t132 * t233 + (t172 * t203 - t174 * t228) * t153) * t148) * t142, 0, 0;];
JaD_rot  = t1;
