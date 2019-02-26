% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPP1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:52
% EndTime: 2019-02-26 21:34:53
% DurationCPUTime: 1.04s
% Computational Cost: add. (6874->125), mult. (6168->272), div. (1114->15), fcn. (7752->9), ass. (0->115)
t179 = qJ(2) + pkin(9);
t177 = cos(t179);
t178 = qJ(4) + pkin(10);
t174 = sin(t178);
t252 = sin(qJ(1));
t213 = t252 * t174;
t176 = cos(t178);
t181 = cos(qJ(1));
t233 = t181 * t176;
t156 = t177 * t233 + t213;
t150 = 0.1e1 / t156 ^ 2;
t175 = sin(t179);
t170 = t175 ^ 2;
t180 = t181 ^ 2;
t237 = t170 * t180;
t218 = t150 * t237;
t146 = 0.1e1 + t218;
t205 = qJD(1) * t252;
t230 = qJD(2) * t181;
t209 = t175 * t230;
t191 = t177 * t205 + t209;
t204 = t252 * qJD(4);
t234 = t181 * t174;
t135 = (-qJD(4) * t177 + qJD(1)) * t234 + (t204 - t191) * t176;
t149 = 0.1e1 / t156;
t247 = t135 * t149 * t150;
t199 = t237 * t247;
t210 = qJD(2) * t175 * t180;
t255 = (-t199 + (-t170 * t181 * t205 + t177 * t210) * t150) / t146 ^ 2;
t235 = t175 * t181;
t152 = t177 * t213 + t233;
t196 = t174 * t204;
t227 = qJD(4) * t181;
t207 = t176 * t227;
t134 = t152 * qJD(1) + t174 * t209 - t177 * t207 - t196;
t212 = t252 * t176;
t155 = t177 * t234 - t212;
t167 = 0.1e1 / t174;
t168 = 0.1e1 / t174 ^ 2;
t171 = 0.1e1 / t175;
t172 = 0.1e1 / t175 ^ 2;
t231 = qJD(2) * t177;
t211 = t172 * t231;
t228 = qJD(4) * t176;
t240 = t167 * t171;
t254 = (t168 * t171 * t228 + t167 * t211) * t155 + t134 * t240;
t236 = t175 * t174;
t142 = atan2(-t152, t236);
t139 = cos(t142);
t138 = sin(t142);
t246 = t138 * t152;
t133 = t139 * t236 - t246;
t130 = 0.1e1 / t133;
t131 = 0.1e1 / t133 ^ 2;
t253 = 0.2e1 * t155;
t147 = t152 ^ 2;
t239 = t168 * t172;
t143 = t147 * t239 + 0.1e1;
t140 = 0.1e1 / t143;
t192 = t174 * t231 + t175 * t228;
t216 = t152 * t239;
t214 = t175 * t252;
t197 = qJD(2) * t214;
t198 = t176 * t205;
t232 = qJD(1) * t181;
t136 = t176 * t204 * t177 - t198 + (t232 * t177 - t197 - t227) * t174;
t219 = t136 * t240;
t122 = (t192 * t216 - t219) * t140;
t189 = -t122 * t152 + t192;
t118 = (-t122 * t236 - t136) * t138 + t189 * t139;
t132 = t130 * t131;
t251 = t118 * t132;
t169 = t167 * t168;
t173 = t171 / t170;
t208 = t172 * t228;
t250 = (t136 * t216 + (-t168 * t173 * t231 - t169 * t208) * t147) / t143 ^ 2;
t249 = t131 * t134;
t248 = t131 * t155;
t245 = t138 * t155;
t244 = t138 * t175;
t243 = t139 * t152;
t242 = t139 * t155;
t241 = t139 * t177;
t238 = t168 * t176;
t229 = qJD(4) * t174;
t148 = t155 ^ 2;
t128 = t131 * t148 + 0.1e1;
t226 = 0.2e1 * (-t134 * t248 - t148 * t251) / t128 ^ 2;
t225 = -0.2e1 * t250;
t224 = 0.2e1 * t255;
t223 = t132 * t253;
t222 = t171 * t250;
t221 = t131 * t245;
t217 = t152 * t240;
t215 = t167 * t172 * t177;
t194 = t152 * t215 + t252;
t129 = t194 * t140;
t206 = t252 - t129;
t203 = t130 * t226;
t202 = t131 * t226;
t201 = t235 * t253;
t200 = t167 * t222;
t154 = t177 * t212 - t234;
t195 = t152 * t238 - t154 * t167;
t193 = t150 * t154 * t181 - t252 * t149;
t144 = 0.1e1 / t146;
t137 = t156 * qJD(1) - t176 * t197 - t177 * t196 - t207;
t126 = 0.1e1 / t128;
t125 = t195 * t171 * t140;
t121 = (-t138 + (t139 * t217 + t138) * t140) * t155;
t120 = -t129 * t243 + (t206 * t244 + t241) * t174;
t119 = t139 * t175 * t176 - t138 * t154 + (-t138 * t236 - t243) * t125;
t117 = t194 * t225 + (t136 * t215 + t232 + (-t168 * t177 * t208 + (-0.2e1 * t173 * t177 ^ 2 - t171) * t167 * qJD(2)) * t152) * t140;
t115 = -0.2e1 * t195 * t222 + (-t195 * t211 + (t136 * t238 - t137 * t167 + (t154 * t238 + (-0.2e1 * t169 * t176 ^ 2 - t167) * t152) * qJD(4)) * t171) * t140;
t1 = [t254 * t140 + t200 * t253, t117, 0, t115, 0, 0; t152 * t203 + (-t136 * t130 + (t118 * t152 + t121 * t134) * t131) * t126 + (t121 * t202 + (0.2e1 * t121 * t251 + (t134 * t140 - t134 - (-t122 * t140 * t217 + t225) * t155) * t131 * t138 + (-(-0.2e1 * t152 * t200 - t122) * t248 + (-(t122 + t219) * t155 + t254 * t152) * t131 * t140) * t139) * t126) * t155, t120 * t155 * t202 + (-(-t117 * t243 + (t122 * t246 - t136 * t139) * t129) * t248 + (t118 * t223 + t249) * t120 + (-t130 * t235 - (-t129 * t244 + t138 * t214 + t241) * t248) * t228) * t126 + (t203 * t235 + ((-t130 * t230 - (t206 * qJD(2) - t122) * t221) * t177 + (t130 * t205 + (t181 * t118 - (-t117 + t232) * t245 - (t206 * t122 - qJD(2)) * t242) * t131) * t175) * t126) * t174, 0 (t119 * t248 - t130 * t156) * t226 + (t119 * t249 + t135 * t130 + (t119 * t223 - t131 * t156) * t118 - (t176 * t231 - t175 * t229 - t115 * t152 - t125 * t136 + (-t125 * t236 - t154) * t122) * t131 * t242 - (-t137 + (-t115 * t174 - t122 * t176) * t175 - t189 * t125) * t221) * t126, 0, 0; t193 * t175 * t224 + (-t193 * t231 + ((qJD(1) * t149 + 0.2e1 * t154 * t247) * t181 + (-t252 * t135 - t137 * t181 + t154 * t205) * t150) * t175) * t144 (t149 * t177 * t181 + t176 * t218) * t224 + (0.2e1 * t176 * t199 + t191 * t149 + ((t135 * t181 - 0.2e1 * t176 * t210) * t177 + (t180 * t229 + 0.2e1 * t181 * t198) * t170) * t150) * t144, 0, t150 * t201 * t255 + (t201 * t247 + (t134 * t235 + (t175 * t205 - t177 * t230) * t155) * t150) * t144, 0, 0;];
JaD_rot  = t1;
