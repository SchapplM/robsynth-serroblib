% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP5
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
% Datum: 2019-02-26 20:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:45:53
% EndTime: 2019-02-26 20:45:54
% DurationCPUTime: 1.09s
% Computational Cost: add. (6874->125), mult. (6168->272), div. (1114->15), fcn. (7752->9), ass. (0->115)
t173 = pkin(9) + qJ(3);
t171 = cos(t173);
t172 = pkin(10) + qJ(5);
t168 = sin(t172);
t246 = sin(qJ(1));
t207 = t246 * t168;
t170 = cos(t172);
t175 = cos(qJ(1));
t227 = t175 * t170;
t150 = t171 * t227 + t207;
t144 = 0.1e1 / t150 ^ 2;
t169 = sin(t173);
t164 = t169 ^ 2;
t174 = t175 ^ 2;
t231 = t164 * t174;
t212 = t144 * t231;
t140 = 0.1e1 + t212;
t199 = qJD(1) * t246;
t224 = qJD(3) * t175;
t203 = t169 * t224;
t185 = t171 * t199 + t203;
t198 = t246 * qJD(5);
t228 = t175 * t168;
t129 = (-qJD(5) * t171 + qJD(1)) * t228 + (t198 - t185) * t170;
t143 = 0.1e1 / t150;
t241 = t129 * t143 * t144;
t193 = t231 * t241;
t204 = qJD(3) * t169 * t174;
t249 = (-t193 + (-t164 * t175 * t199 + t171 * t204) * t144) / t140 ^ 2;
t229 = t169 * t175;
t146 = t171 * t207 + t227;
t190 = t168 * t198;
t221 = qJD(5) * t175;
t201 = t170 * t221;
t128 = t146 * qJD(1) + t168 * t203 - t171 * t201 - t190;
t206 = t246 * t170;
t149 = t171 * t228 - t206;
t161 = 0.1e1 / t168;
t162 = 0.1e1 / t168 ^ 2;
t165 = 0.1e1 / t169;
t166 = 0.1e1 / t169 ^ 2;
t225 = qJD(3) * t171;
t205 = t166 * t225;
t222 = qJD(5) * t170;
t234 = t161 * t165;
t248 = (t162 * t165 * t222 + t161 * t205) * t149 + t128 * t234;
t230 = t169 * t168;
t136 = atan2(-t146, t230);
t133 = cos(t136);
t132 = sin(t136);
t240 = t132 * t146;
t127 = t133 * t230 - t240;
t124 = 0.1e1 / t127;
t125 = 0.1e1 / t127 ^ 2;
t247 = 0.2e1 * t149;
t141 = t146 ^ 2;
t233 = t162 * t166;
t137 = t141 * t233 + 0.1e1;
t134 = 0.1e1 / t137;
t186 = t168 * t225 + t169 * t222;
t210 = t146 * t233;
t208 = t169 * t246;
t191 = qJD(3) * t208;
t192 = t170 * t199;
t226 = qJD(1) * t175;
t130 = t170 * t198 * t171 - t192 + (t226 * t171 - t191 - t221) * t168;
t213 = t130 * t234;
t116 = (t186 * t210 - t213) * t134;
t183 = -t116 * t146 + t186;
t112 = (-t116 * t230 - t130) * t132 + t183 * t133;
t126 = t124 * t125;
t245 = t112 * t126;
t163 = t161 * t162;
t167 = t165 / t164;
t202 = t166 * t222;
t244 = (t130 * t210 + (-t162 * t167 * t225 - t163 * t202) * t141) / t137 ^ 2;
t243 = t125 * t149;
t242 = t128 * t125;
t239 = t132 * t149;
t238 = t132 * t169;
t237 = t133 * t146;
t236 = t133 * t149;
t235 = t133 * t171;
t232 = t162 * t170;
t223 = qJD(5) * t168;
t142 = t149 ^ 2;
t122 = t125 * t142 + 0.1e1;
t220 = 0.2e1 * (-t142 * t245 - t149 * t242) / t122 ^ 2;
t219 = -0.2e1 * t244;
t218 = 0.2e1 * t249;
t217 = t126 * t247;
t216 = t165 * t244;
t215 = t125 * t239;
t211 = t146 * t234;
t209 = t161 * t166 * t171;
t188 = t146 * t209 + t246;
t123 = t188 * t134;
t200 = t246 - t123;
t197 = t124 * t220;
t196 = t125 * t220;
t195 = t229 * t247;
t194 = t161 * t216;
t148 = t171 * t206 - t228;
t189 = t146 * t232 - t148 * t161;
t187 = t144 * t148 * t175 - t246 * t143;
t138 = 0.1e1 / t140;
t131 = t150 * qJD(1) - t170 * t191 - t171 * t190 - t201;
t120 = 0.1e1 / t122;
t119 = t189 * t165 * t134;
t115 = (-t132 + (t133 * t211 + t132) * t134) * t149;
t114 = -t123 * t237 + (t200 * t238 + t235) * t168;
t113 = t133 * t169 * t170 - t132 * t148 + (-t132 * t230 - t237) * t119;
t111 = t188 * t219 + (t130 * t209 + t226 + (-t162 * t171 * t202 + (-0.2e1 * t167 * t171 ^ 2 - t165) * t161 * qJD(3)) * t146) * t134;
t109 = -0.2e1 * t189 * t216 + (-t189 * t205 + (t130 * t232 - t131 * t161 + (t148 * t232 + (-0.2e1 * t163 * t170 ^ 2 - t161) * t146) * qJD(5)) * t165) * t134;
t1 = [t248 * t134 + t194 * t247, 0, t111, 0, t109, 0; t146 * t197 + (-t130 * t124 + (t112 * t146 + t115 * t128) * t125) * t120 + (t115 * t196 + (0.2e1 * t115 * t245 + (t128 * t134 - t128 - (-t116 * t134 * t211 + t219) * t149) * t125 * t132 + (-(-0.2e1 * t146 * t194 - t116) * t243 + (-(t116 + t213) * t149 + t248 * t146) * t125 * t134) * t133) * t120) * t149, 0, t114 * t149 * t196 + (-(-t111 * t237 + (t116 * t240 - t130 * t133) * t123) * t243 + (t112 * t217 + t242) * t114 + (-t124 * t229 - (-t123 * t238 + t132 * t208 + t235) * t243) * t222) * t120 + (t197 * t229 + ((-t124 * t224 - (t200 * qJD(3) - t116) * t215) * t171 + (t124 * t199 + (t175 * t112 - (-t111 + t226) * t239 - (t200 * t116 - qJD(3)) * t236) * t125) * t169) * t120) * t168, 0 (t113 * t243 - t124 * t150) * t220 + (t113 * t242 + t129 * t124 + (t113 * t217 - t150 * t125) * t112 - (t170 * t225 - t169 * t223 - t109 * t146 - t119 * t130 + (-t119 * t230 - t148) * t116) * t125 * t236 - (-t131 + (-t109 * t168 - t116 * t170) * t169 - t183 * t119) * t215) * t120, 0; t187 * t169 * t218 + (-t187 * t225 + ((qJD(1) * t143 + 0.2e1 * t148 * t241) * t175 + (-t246 * t129 - t131 * t175 + t148 * t199) * t144) * t169) * t138, 0 (t143 * t171 * t175 + t170 * t212) * t218 + (0.2e1 * t170 * t193 + t185 * t143 + ((t129 * t175 - 0.2e1 * t170 * t204) * t171 + (t174 * t223 + 0.2e1 * t175 * t192) * t164) * t144) * t138, 0, t144 * t195 * t249 + (t195 * t241 + (t128 * t229 + (t169 * t199 - t171 * t224) * t149) * t144) * t138, 0;];
JaD_rot  = t1;
