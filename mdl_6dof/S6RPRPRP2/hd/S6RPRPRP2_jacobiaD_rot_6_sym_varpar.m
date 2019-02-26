% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:44:15
% EndTime: 2019-02-26 20:44:16
% DurationCPUTime: 1.01s
% Computational Cost: add. (6087->126), mult. (6168->275), div. (1114->15), fcn. (7752->9), ass. (0->116)
t171 = qJ(3) + pkin(10);
t169 = cos(t171);
t175 = sin(qJ(5));
t223 = qJ(1) + pkin(9);
t205 = sin(t223);
t198 = t205 * t175;
t170 = cos(t223);
t176 = cos(qJ(5));
t231 = t170 * t176;
t152 = t169 * t231 + t198;
t146 = 0.1e1 / t152 ^ 2;
t168 = sin(t171);
t163 = t168 ^ 2;
t167 = t170 ^ 2;
t238 = t163 * t167;
t214 = t146 * t238;
t136 = 0.1e1 + t214;
t195 = qJD(1) * t205;
t191 = t169 * t195;
t194 = t205 * qJD(5);
t226 = qJD(3) * t176;
t208 = t168 * t226;
t131 = (-t191 + t194) * t176 + (-t208 + (-qJD(5) * t169 + qJD(1)) * t175) * t170;
t145 = 0.1e1 / t152;
t245 = t131 * t145 * t146;
t200 = t238 * t245;
t228 = qJD(3) * t169;
t252 = (-t200 + (-t163 * t170 * t195 + t167 * t168 * t228) * t146) / t136 ^ 2;
t234 = t168 * t170;
t148 = t169 * t198 + t231;
t190 = t175 * t194;
t224 = qJD(5) * t176;
t206 = t170 * t224;
t227 = qJD(3) * t170;
t209 = t168 * t227;
t130 = t148 * qJD(1) - t169 * t206 + t175 * t209 - t190;
t197 = t205 * t176;
t232 = t170 * t175;
t151 = t169 * t232 - t197;
t164 = 0.1e1 / t168;
t172 = 0.1e1 / t175;
t173 = 0.1e1 / t175 ^ 2;
t207 = t173 * t224;
t165 = 0.1e1 / t168 ^ 2;
t210 = t165 * t228;
t237 = t164 * t172;
t251 = (t164 * t207 + t172 * t210) * t151 + t130 * t237;
t233 = t168 * t175;
t141 = atan2(-t148, t233);
t138 = cos(t141);
t137 = sin(t141);
t244 = t137 * t148;
t129 = t138 * t233 - t244;
t126 = 0.1e1 / t129;
t127 = 0.1e1 / t129 ^ 2;
t250 = 0.2e1 * t151;
t143 = t148 ^ 2;
t235 = t165 * t173;
t142 = t143 * t235 + 0.1e1;
t139 = 0.1e1 / t142;
t188 = t168 * t224 + t175 * t228;
t212 = t148 * t235;
t199 = t168 * t205;
t192 = qJD(3) * t199;
t193 = t176 * t195;
t225 = qJD(5) * t175;
t229 = qJD(1) * t170;
t132 = -t175 * t192 - t170 * t225 - t193 + (t175 * t229 + t176 * t194) * t169;
t215 = t132 * t237;
t118 = (t188 * t212 - t215) * t139;
t184 = -t118 * t148 + t188;
t114 = (-t118 * t233 - t132) * t137 + t184 * t138;
t128 = t126 * t127;
t249 = t114 * t128;
t166 = t164 / t163;
t174 = t172 * t173;
t248 = (t132 * t212 + (-t165 * t174 * t224 - t166 * t173 * t228) * t143) / t142 ^ 2;
t247 = t127 * t151;
t246 = t130 * t127;
t243 = t137 * t151;
t242 = t137 * t168;
t241 = t138 * t148;
t240 = t138 * t151;
t239 = t138 * t169;
t236 = t165 * t169;
t230 = t173 * t176;
t144 = t151 ^ 2;
t124 = t127 * t144 + 0.1e1;
t222 = 0.2e1 * (-t144 * t249 - t151 * t246) / t124 ^ 2;
t221 = 0.2e1 * t252;
t220 = -0.2e1 * t248;
t219 = t128 * t250;
t218 = t164 * t248;
t217 = t127 * t243;
t213 = t148 * t237;
t211 = t172 * t236;
t204 = t126 * t222;
t203 = t127 * t222;
t202 = t234 * t250;
t201 = t172 * t218;
t187 = t148 * t211 + t205;
t125 = t187 * t139;
t196 = t205 - t125;
t150 = t169 * t197 - t232;
t189 = t148 * t230 - t150 * t172;
t186 = t146 * t150 * t170 - t205 * t145;
t134 = 0.1e1 / t136;
t133 = t152 * qJD(1) - t169 * t190 - t176 * t192 - t206;
t122 = 0.1e1 / t124;
t121 = t189 * t164 * t139;
t117 = (-t137 + (t138 * t213 + t137) * t139) * t151;
t116 = -t125 * t241 + (t196 * t242 + t239) * t175;
t115 = t138 * t168 * t176 - t137 * t150 + (-t137 * t233 - t241) * t121;
t113 = t187 * t220 + (t132 * t211 + t229 + (-t207 * t236 + (-0.2e1 * t166 * t169 ^ 2 - t164) * t172 * qJD(3)) * t148) * t139;
t111 = -0.2e1 * t189 * t218 + (-t189 * t210 + (t132 * t230 - t133 * t172 + (t150 * t230 + (-0.2e1 * t174 * t176 ^ 2 - t172) * t148) * qJD(5)) * t164) * t139;
t1 = [t251 * t139 + t201 * t250, 0, t113, 0, t111, 0; t148 * t204 + (-t132 * t126 + (t114 * t148 + t117 * t130) * t127) * t122 + (t117 * t203 + (0.2e1 * t117 * t249 + (t130 * t139 - t130 - (-t118 * t139 * t213 + t220) * t151) * t127 * t137 + (-(-0.2e1 * t148 * t201 - t118) * t247 + (-(t118 + t215) * t151 + t251 * t148) * t127 * t139) * t138) * t122) * t151, 0, t116 * t151 * t203 + (-(-t113 * t241 + (t118 * t244 - t132 * t138) * t125) * t247 + (t114 * t219 + t246) * t116 + (-t126 * t234 - (-t125 * t242 + t137 * t199 + t239) * t247) * t224) * t122 + (t204 * t234 + ((-t126 * t227 - (t196 * qJD(3) - t118) * t217) * t169 + (t126 * t195 + (t170 * t114 - (-t113 + t229) * t243 - (t196 * t118 - qJD(3)) * t240) * t127) * t168) * t122) * t175, 0 (t115 * t247 - t126 * t152) * t222 + (t115 * t246 + t131 * t126 + (t115 * t219 - t152 * t127) * t114 - (t169 * t226 - t168 * t225 - t111 * t148 - t121 * t132 + (-t121 * t233 - t150) * t118) * t127 * t240 - (-t133 + (-t111 * t175 - t118 * t176) * t168 - t184 * t121) * t217) * t122, 0; t186 * t168 * t221 + (-t186 * t228 + ((qJD(1) * t145 + 0.2e1 * t150 * t245) * t170 + (-t205 * t131 - t133 * t170 + t150 * t195) * t146) * t168) * t134, 0 (t145 * t169 * t170 + t176 * t214) * t221 + (0.2e1 * t176 * t200 + (t191 + t209) * t145 + ((t131 * t170 - 0.2e1 * t167 * t208) * t169 + (t167 * t225 + 0.2e1 * t170 * t193) * t163) * t146) * t134, 0, t146 * t202 * t252 + (t202 * t245 + (t130 * t234 + (t168 * t195 - t169 * t227) * t151) * t146) * t134, 0;];
JaD_rot  = t1;
