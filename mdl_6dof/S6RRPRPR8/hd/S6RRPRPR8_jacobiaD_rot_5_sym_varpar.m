% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR8_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:41:47
% EndTime: 2019-02-26 21:41:48
% DurationCPUTime: 0.99s
% Computational Cost: add. (4821->123), mult. (6168->268), div. (1114->15), fcn. (7752->9), ass. (0->114)
t156 = pkin(10) + qJ(4);
t155 = cos(t156);
t238 = 0.2e1 * t155;
t154 = sin(t156);
t163 = cos(qJ(2));
t164 = cos(qJ(1));
t214 = t164 * t155;
t234 = sin(qJ(1));
t140 = t234 * t154 + t163 * t214;
t134 = 0.1e1 / t140 ^ 2;
t162 = sin(qJ(2));
t157 = t162 ^ 2;
t161 = t164 ^ 2;
t219 = t157 * t161;
t199 = t134 * t219;
t130 = 0.1e1 + t199;
t188 = qJD(1) * t234;
t212 = qJD(2) * t163;
t172 = t157 * t164 * t188 - t161 * t162 * t212;
t211 = qJD(2) * t164;
t192 = t162 * t211;
t175 = t163 * t188 + t192;
t187 = t234 * qJD(4);
t215 = t164 * t154;
t119 = (-qJD(4) * t163 + qJD(1)) * t215 + (t187 - t175) * t155;
t133 = 0.1e1 / t140;
t229 = t119 * t133 * t134;
t182 = t219 * t229;
t237 = (-t172 * t134 - t182) / t130 ^ 2;
t217 = t162 * t164;
t194 = t234 * t163;
t136 = t154 * t194 + t214;
t180 = t154 * t187;
t208 = qJD(4) * t164;
t190 = t155 * t208;
t118 = t136 * qJD(1) + t154 * t192 - t163 * t190 - t180;
t139 = -t234 * t155 + t163 * t215;
t151 = 0.1e1 / t154;
t152 = 0.1e1 / t154 ^ 2;
t158 = 0.1e1 / t162;
t159 = 0.1e1 / t162 ^ 2;
t193 = t159 * t212;
t210 = qJD(4) * t155;
t222 = t151 * t158;
t236 = (t152 * t158 * t210 + t151 * t193) * t139 + t118 * t222;
t218 = t162 * t154;
t126 = atan2(-t136, t218);
t123 = cos(t126);
t122 = sin(t126);
t228 = t122 * t136;
t117 = t123 * t218 - t228;
t114 = 0.1e1 / t117;
t115 = 0.1e1 / t117 ^ 2;
t235 = 0.2e1 * t139;
t131 = t136 ^ 2;
t220 = t152 * t159;
t127 = t131 * t220 + 0.1e1;
t124 = 0.1e1 / t127;
t209 = qJD(4) * t162;
t176 = t154 * t212 + t155 * t209;
t197 = t136 * t220;
t195 = t234 * t162;
t181 = qJD(2) * t195;
t213 = qJD(1) * t164;
t120 = (t187 * t163 - t188) * t155 + (t213 * t163 - t181 - t208) * t154;
t200 = t120 * t222;
t106 = (t176 * t197 - t200) * t124;
t173 = -t106 * t136 + t176;
t102 = (-t106 * t218 - t120) * t122 + t173 * t123;
t116 = t114 * t115;
t233 = t102 * t116;
t153 = t151 * t152;
t160 = t158 / t157;
t191 = t159 * t210;
t232 = (t120 * t197 + (-t152 * t160 * t212 - t153 * t191) * t131) / t127 ^ 2;
t231 = t115 * t139;
t230 = t118 * t115;
t227 = t122 * t139;
t226 = t122 * t162;
t225 = t123 * t136;
t224 = t123 * t139;
t223 = t123 * t163;
t221 = t152 * t155;
t216 = t163 * t164;
t132 = t139 ^ 2;
t112 = t115 * t132 + 0.1e1;
t207 = 0.2e1 * (-t132 * t233 - t139 * t230) / t112 ^ 2;
t206 = -0.2e1 * t232;
t205 = 0.2e1 * t237;
t204 = t116 * t235;
t203 = t158 * t232;
t202 = t115 * t227;
t198 = t136 * t222;
t196 = t151 * t159 * t163;
t178 = t136 * t196 + t234;
t113 = t178 * t124;
t189 = t234 - t113;
t186 = t114 * t207;
t185 = t115 * t207;
t184 = t217 * t235;
t183 = t151 * t203;
t138 = t155 * t194 - t215;
t179 = t136 * t221 - t138 * t151;
t177 = t134 * t138 * t164 - t234 * t133;
t128 = 0.1e1 / t130;
t121 = t140 * qJD(1) - t155 * t181 - t163 * t180 - t190;
t110 = 0.1e1 / t112;
t109 = t179 * t158 * t124;
t105 = (-t122 + (t123 * t198 + t122) * t124) * t139;
t104 = -t113 * t225 + (t189 * t226 + t223) * t154;
t103 = t123 * t155 * t162 - t122 * t138 + (-t122 * t218 - t225) * t109;
t101 = t178 * t206 + (t120 * t196 + t213 + (-t152 * t163 * t191 + (-0.2e1 * t160 * t163 ^ 2 - t158) * t151 * qJD(2)) * t136) * t124;
t99 = -0.2e1 * t179 * t203 + (-t179 * t193 + (t120 * t221 - t121 * t151 + (t138 * t221 + (-0.2e1 * t153 * t155 ^ 2 - t151) * t136) * qJD(4)) * t158) * t124;
t1 = [t236 * t124 + t183 * t235, t101, 0, t99, 0, 0; t136 * t186 + (-t120 * t114 + (t102 * t136 + t105 * t118) * t115) * t110 + (t105 * t185 + (0.2e1 * t105 * t233 + (t118 * t124 - t118 - (-t106 * t124 * t198 + t206) * t139) * t115 * t122 + (-(-0.2e1 * t136 * t183 - t106) * t231 + (-(t106 + t200) * t139 + t236 * t136) * t115 * t124) * t123) * t110) * t139, t104 * t139 * t185 + (-(-t101 * t225 + (t106 * t228 - t120 * t123) * t113) * t231 + (t102 * t204 + t230) * t104 + (-t114 * t217 - (-t113 * t226 + t122 * t195 + t223) * t231) * t210) * t110 + (t186 * t217 + ((-t114 * t211 - (t189 * qJD(2) - t106) * t202) * t163 + (t114 * t188 + (t164 * t102 - (-t101 + t213) * t227 - (t189 * t106 - qJD(2)) * t224) * t115) * t162) * t110) * t154, 0 (t103 * t231 - t114 * t140) * t207 + (t103 * t230 + t119 * t114 + (t103 * t204 - t115 * t140) * t102 - (t155 * t212 - t154 * t209 - t109 * t120 - t136 * t99 + (-t109 * t218 - t138) * t106) * t115 * t224 - (-t121 + (-t106 * t155 - t154 * t99) * t162 - t173 * t109) * t202) * t110, 0, 0; t177 * t162 * t205 + (-t177 * t212 + ((qJD(1) * t133 + 0.2e1 * t138 * t229) * t164 + (-t234 * t119 - t121 * t164 + t138 * t188) * t134) * t162) * t128 (t133 * t216 + t155 * t199) * t205 + (t182 * t238 + t175 * t133 + (qJD(4) * t154 * t219 + t119 * t216 + t172 * t238) * t134) * t128, 0, t134 * t184 * t237 + (t184 * t229 + (t118 * t217 + (t162 * t188 - t163 * t211) * t139) * t134) * t128, 0, 0;];
JaD_rot  = t1;
