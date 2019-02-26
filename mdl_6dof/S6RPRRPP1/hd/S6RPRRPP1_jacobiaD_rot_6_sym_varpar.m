% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP1
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
% Datum: 2019-02-26 20:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:56:17
% EndTime: 2019-02-26 20:56:18
% DurationCPUTime: 1.07s
% Computational Cost: add. (7031->125), mult. (6168->273), div. (1114->15), fcn. (7752->9), ass. (0->115)
t164 = qJ(4) + pkin(10);
t162 = cos(t164);
t161 = sin(t164);
t217 = qJ(1) + pkin(9);
t199 = sin(t217);
t192 = t199 * t161;
t163 = cos(t217);
t170 = cos(qJ(3));
t225 = t163 * t170;
t146 = t162 * t225 + t192;
t140 = 0.1e1 / t146 ^ 2;
t160 = t163 ^ 2;
t169 = sin(qJ(3));
t165 = t169 ^ 2;
t228 = t160 * t165;
t208 = t140 * t228;
t135 = 0.1e1 + t208;
t189 = qJD(1) * t199;
t222 = qJD(3) * t163;
t204 = t169 * t222;
t178 = t170 * t189 + t204;
t188 = t199 * qJD(4);
t227 = t163 * t161;
t125 = (-qJD(4) * t170 + qJD(1)) * t227 + (t188 - t178) * t162;
t139 = 0.1e1 / t146;
t238 = t125 * t139 * t140;
t194 = t228 * t238;
t221 = qJD(3) * t170;
t202 = t169 * t221;
t245 = (-t194 + (-t163 * t165 * t189 + t160 * t202) * t140) / t135 ^ 2;
t226 = t163 * t169;
t142 = t163 * t162 + t170 * t192;
t185 = t161 * t188;
t219 = qJD(4) * t163;
t200 = t162 * t219;
t124 = t142 * qJD(1) + t161 * t204 - t170 * t200 - t185;
t191 = t199 * t162;
t145 = t161 * t225 - t191;
t157 = 0.1e1 / t161;
t158 = 0.1e1 / t161 ^ 2;
t166 = 0.1e1 / t169;
t167 = 0.1e1 / t169 ^ 2;
t203 = t167 * t221;
t220 = qJD(4) * t162;
t231 = t157 * t166;
t244 = (t158 * t166 * t220 + t157 * t203) * t145 + t124 * t231;
t224 = t169 * t161;
t134 = atan2(-t142, t224);
t129 = cos(t134);
t128 = sin(t134);
t237 = t128 * t142;
t123 = t129 * t224 - t237;
t120 = 0.1e1 / t123;
t121 = 0.1e1 / t123 ^ 2;
t243 = 0.2e1 * t145;
t137 = t142 ^ 2;
t229 = t158 * t167;
t136 = t137 * t229 + 0.1e1;
t132 = 0.1e1 / t136;
t218 = qJD(4) * t169;
t183 = t161 * t221 + t162 * t218;
t206 = t142 * t229;
t193 = t169 * t199;
t186 = qJD(3) * t193;
t187 = t162 * t189;
t223 = qJD(1) * t163;
t126 = t162 * t188 * t170 - t187 + (t223 * t170 - t186 - t219) * t161;
t209 = t126 * t231;
t112 = (t183 * t206 - t209) * t132;
t179 = -t112 * t142 + t183;
t108 = (-t112 * t224 - t126) * t128 + t179 * t129;
t122 = t120 * t121;
t242 = t108 * t122;
t159 = t157 * t158;
t168 = t166 / t165;
t201 = t167 * t220;
t241 = (t126 * t206 + (-t158 * t168 * t221 - t159 * t201) * t137) / t136 ^ 2;
t240 = t121 * t145;
t239 = t124 * t121;
t236 = t128 * t145;
t235 = t128 * t169;
t234 = t129 * t142;
t233 = t129 * t145;
t232 = t129 * t170;
t230 = t158 * t162;
t138 = t145 ^ 2;
t118 = t121 * t138 + 0.1e1;
t216 = 0.2e1 * (-t138 * t242 - t145 * t239) / t118 ^ 2;
t215 = 0.2e1 * t245;
t214 = -0.2e1 * t241;
t213 = t122 * t243;
t212 = t166 * t241;
t211 = t121 * t236;
t207 = t142 * t231;
t205 = t157 * t167 * t170;
t198 = t120 * t216;
t197 = t121 * t216;
t196 = t226 * t243;
t195 = t157 * t212;
t182 = t142 * t205 + t199;
t119 = t182 * t132;
t190 = t199 - t119;
t144 = t170 * t191 - t227;
t184 = t142 * t230 - t144 * t157;
t181 = t140 * t144 * t163 - t199 * t139;
t130 = 0.1e1 / t135;
t127 = t146 * qJD(1) - t162 * t186 - t170 * t185 - t200;
t116 = 0.1e1 / t118;
t115 = t184 * t166 * t132;
t111 = (-t128 + (t129 * t207 + t128) * t132) * t145;
t110 = -t119 * t234 + (t190 * t235 + t232) * t161;
t109 = t129 * t162 * t169 - t128 * t144 + (-t128 * t224 - t234) * t115;
t107 = t182 * t214 + (t126 * t205 + t223 + (-t158 * t170 * t201 + (-0.2e1 * t168 * t170 ^ 2 - t166) * t157 * qJD(3)) * t142) * t132;
t105 = -0.2e1 * t184 * t212 + (-t184 * t203 + (t126 * t230 - t127 * t157 + (t144 * t230 + (-0.2e1 * t159 * t162 ^ 2 - t157) * t142) * qJD(4)) * t166) * t132;
t1 = [t244 * t132 + t195 * t243, 0, t107, t105, 0, 0; t142 * t198 + (-t126 * t120 + (t108 * t142 + t111 * t124) * t121) * t116 + (t111 * t197 + (0.2e1 * t111 * t242 + (t124 * t132 - t124 - (-t112 * t132 * t207 + t214) * t145) * t121 * t128 + (-(-0.2e1 * t142 * t195 - t112) * t240 + (-(t112 + t209) * t145 + t244 * t142) * t121 * t132) * t129) * t116) * t145, 0, t110 * t145 * t197 + (-(-t107 * t234 + (t112 * t237 - t126 * t129) * t119) * t240 + (t108 * t213 + t239) * t110 + (-t120 * t226 - (-t119 * t235 + t128 * t193 + t232) * t240) * t220) * t116 + (t198 * t226 + ((-t120 * t222 - (t190 * qJD(3) - t112) * t211) * t170 + (t120 * t189 + (t163 * t108 - (-t107 + t223) * t236 - (t190 * t112 - qJD(3)) * t233) * t121) * t169) * t116) * t161 (t109 * t240 - t120 * t146) * t216 + (t109 * t239 + t125 * t120 + (t109 * t213 - t146 * t121) * t108 - (t162 * t221 - t161 * t218 - t105 * t142 - t115 * t126 + (-t115 * t224 - t144) * t112) * t121 * t233 - (-t127 + (-t105 * t161 - t112 * t162) * t169 - t179 * t115) * t211) * t116, 0, 0; t181 * t169 * t215 + (-t181 * t221 + ((qJD(1) * t139 + 0.2e1 * t144 * t238) * t163 + (-t199 * t125 - t127 * t163 + t144 * t189) * t140) * t169) * t130, 0 (t139 * t225 + t162 * t208) * t215 + (0.2e1 * t162 * t194 + t178 * t139 + ((t125 * t170 + 0.2e1 * t165 * t187) * t163 + (qJD(4) * t161 * t165 - 0.2e1 * t162 * t202) * t160) * t140) * t130, t140 * t196 * t245 + (t196 * t238 + (t124 * t226 + (-t163 * t221 + t169 * t189) * t145) * t140) * t130, 0, 0;];
JaD_rot  = t1;
