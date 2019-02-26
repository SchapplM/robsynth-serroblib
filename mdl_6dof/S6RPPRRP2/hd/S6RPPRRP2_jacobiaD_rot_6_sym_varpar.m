% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRP2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:30:41
% EndTime: 2019-02-26 20:30:42
% DurationCPUTime: 1.01s
% Computational Cost: add. (6087->126), mult. (6168->275), div. (1114->15), fcn. (7752->9), ass. (0->116)
t168 = pkin(10) + qJ(4);
t166 = cos(t168);
t172 = sin(qJ(5));
t220 = qJ(1) + pkin(9);
t202 = sin(t220);
t195 = t202 * t172;
t167 = cos(t220);
t173 = cos(qJ(5));
t228 = t167 * t173;
t149 = t166 * t228 + t195;
t143 = 0.1e1 / t149 ^ 2;
t165 = sin(t168);
t160 = t165 ^ 2;
t164 = t167 ^ 2;
t235 = t160 * t164;
t211 = t143 * t235;
t133 = 0.1e1 + t211;
t192 = qJD(1) * t202;
t188 = t166 * t192;
t191 = t202 * qJD(5);
t223 = qJD(4) * t173;
t205 = t165 * t223;
t128 = (-t188 + t191) * t173 + (-t205 + (-qJD(5) * t166 + qJD(1)) * t172) * t167;
t142 = 0.1e1 / t149;
t242 = t128 * t142 * t143;
t197 = t235 * t242;
t225 = qJD(4) * t166;
t249 = (-t197 + (-t160 * t167 * t192 + t164 * t165 * t225) * t143) / t133 ^ 2;
t231 = t165 * t167;
t145 = t166 * t195 + t228;
t187 = t172 * t191;
t221 = qJD(5) * t173;
t203 = t167 * t221;
t224 = qJD(4) * t167;
t206 = t165 * t224;
t127 = t145 * qJD(1) - t166 * t203 + t172 * t206 - t187;
t194 = t202 * t173;
t229 = t167 * t172;
t148 = t166 * t229 - t194;
t161 = 0.1e1 / t165;
t169 = 0.1e1 / t172;
t170 = 0.1e1 / t172 ^ 2;
t204 = t170 * t221;
t162 = 0.1e1 / t165 ^ 2;
t207 = t162 * t225;
t234 = t161 * t169;
t248 = (t161 * t204 + t169 * t207) * t148 + t127 * t234;
t230 = t165 * t172;
t138 = atan2(-t145, t230);
t135 = cos(t138);
t134 = sin(t138);
t241 = t134 * t145;
t126 = t135 * t230 - t241;
t123 = 0.1e1 / t126;
t124 = 0.1e1 / t126 ^ 2;
t247 = 0.2e1 * t148;
t140 = t145 ^ 2;
t232 = t162 * t170;
t139 = t140 * t232 + 0.1e1;
t136 = 0.1e1 / t139;
t185 = t165 * t221 + t172 * t225;
t209 = t145 * t232;
t196 = t165 * t202;
t189 = qJD(4) * t196;
t190 = t173 * t192;
t222 = qJD(5) * t172;
t226 = qJD(1) * t167;
t129 = -t172 * t189 - t167 * t222 - t190 + (t172 * t226 + t173 * t191) * t166;
t212 = t129 * t234;
t115 = (t185 * t209 - t212) * t136;
t181 = -t115 * t145 + t185;
t111 = (-t115 * t230 - t129) * t134 + t181 * t135;
t125 = t123 * t124;
t246 = t111 * t125;
t163 = t161 / t160;
t171 = t169 * t170;
t245 = (t129 * t209 + (-t162 * t171 * t221 - t163 * t170 * t225) * t140) / t139 ^ 2;
t244 = t124 * t148;
t243 = t127 * t124;
t240 = t134 * t148;
t239 = t134 * t165;
t238 = t135 * t145;
t237 = t135 * t148;
t236 = t135 * t166;
t233 = t162 * t166;
t227 = t170 * t173;
t141 = t148 ^ 2;
t121 = t124 * t141 + 0.1e1;
t219 = 0.2e1 * (-t141 * t246 - t148 * t243) / t121 ^ 2;
t218 = 0.2e1 * t249;
t217 = -0.2e1 * t245;
t216 = t125 * t247;
t215 = t161 * t245;
t214 = t124 * t240;
t210 = t145 * t234;
t208 = t169 * t233;
t201 = t123 * t219;
t200 = t124 * t219;
t199 = t231 * t247;
t198 = t169 * t215;
t184 = t145 * t208 + t202;
t122 = t184 * t136;
t193 = t202 - t122;
t147 = t166 * t194 - t229;
t186 = t145 * t227 - t147 * t169;
t183 = t143 * t147 * t167 - t202 * t142;
t131 = 0.1e1 / t133;
t130 = t149 * qJD(1) - t166 * t187 - t173 * t189 - t203;
t119 = 0.1e1 / t121;
t118 = t186 * t161 * t136;
t114 = (-t134 + (t135 * t210 + t134) * t136) * t148;
t113 = -t122 * t238 + (t193 * t239 + t236) * t172;
t112 = t135 * t165 * t173 - t134 * t147 + (-t134 * t230 - t238) * t118;
t110 = t184 * t217 + (t129 * t208 + t226 + (-t204 * t233 + (-0.2e1 * t163 * t166 ^ 2 - t161) * t169 * qJD(4)) * t145) * t136;
t108 = -0.2e1 * t186 * t215 + (-t186 * t207 + (t129 * t227 - t130 * t169 + (t147 * t227 + (-0.2e1 * t171 * t173 ^ 2 - t169) * t145) * qJD(5)) * t161) * t136;
t1 = [t248 * t136 + t198 * t247, 0, 0, t110, t108, 0; t145 * t201 + (-t129 * t123 + (t111 * t145 + t114 * t127) * t124) * t119 + (t114 * t200 + (0.2e1 * t114 * t246 + (t127 * t136 - t127 - (-t115 * t136 * t210 + t217) * t148) * t124 * t134 + (-(-0.2e1 * t145 * t198 - t115) * t244 + (-(t115 + t212) * t148 + t248 * t145) * t124 * t136) * t135) * t119) * t148, 0, 0, t113 * t148 * t200 + (-(-t110 * t238 + (t115 * t241 - t129 * t135) * t122) * t244 + (t111 * t216 + t243) * t113 + (-t123 * t231 - (-t122 * t239 + t134 * t196 + t236) * t244) * t221) * t119 + (t201 * t231 + ((-t123 * t224 - (t193 * qJD(4) - t115) * t214) * t166 + (t123 * t192 + (t167 * t111 - (-t110 + t226) * t240 - (t193 * t115 - qJD(4)) * t237) * t124) * t165) * t119) * t172 (t112 * t244 - t123 * t149) * t219 + (t112 * t243 + t128 * t123 + (t112 * t216 - t149 * t124) * t111 - (t166 * t223 - t165 * t222 - t108 * t145 - t118 * t129 + (-t118 * t230 - t147) * t115) * t124 * t237 - (-t130 + (-t108 * t172 - t115 * t173) * t165 - t181 * t118) * t214) * t119, 0; t183 * t165 * t218 + (-t183 * t225 + ((qJD(1) * t142 + 0.2e1 * t147 * t242) * t167 + (-t202 * t128 - t130 * t167 + t147 * t192) * t143) * t165) * t131, 0, 0 (t142 * t166 * t167 + t173 * t211) * t218 + (0.2e1 * t173 * t197 + (t188 + t206) * t142 + ((t128 * t167 - 0.2e1 * t164 * t205) * t166 + (t164 * t222 + 0.2e1 * t167 * t190) * t160) * t143) * t131, t143 * t199 * t249 + (t199 * t242 + (t127 * t231 + (t165 * t192 - t166 * t224) * t148) * t143) * t131, 0;];
JaD_rot  = t1;
