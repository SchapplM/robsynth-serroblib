% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:41:59
% EndTime: 2019-02-26 20:41:59
% DurationCPUTime: 0.71s
% Computational Cost: add. (2429->94), mult. (2519->204), div. (480->12), fcn. (2968->9), ass. (0->95)
t144 = qJ(3) + pkin(9);
t140 = sin(t144);
t135 = 0.1e1 / t140 ^ 2;
t142 = cos(t144);
t138 = t142 ^ 2;
t190 = t135 * t138;
t148 = cos(qJ(1));
t209 = 0.2e1 * t148;
t208 = t142 * t190;
t184 = t148 * t142;
t129 = atan2(-t184, t140);
t127 = sin(t129);
t128 = cos(t129);
t113 = -t127 * t184 + t128 * t140;
t110 = 0.1e1 / t113;
t143 = pkin(10) + qJ(6);
t139 = sin(t143);
t141 = cos(t143);
t147 = sin(qJ(1));
t186 = t147 * t141;
t124 = t139 * t148 + t140 * t186;
t120 = 0.1e1 / t124;
t134 = 0.1e1 / t140;
t111 = 0.1e1 / t113 ^ 2;
t121 = 0.1e1 / t124 ^ 2;
t146 = t148 ^ 2;
t132 = t146 * t190 + 0.1e1;
t130 = 0.1e1 / t132;
t207 = t130 - 0.1e1;
t145 = t147 ^ 2;
t189 = t138 * t145;
t109 = t111 * t189 + 0.1e1;
t182 = qJD(1) * t148;
t164 = t138 * t147 * t182;
t181 = qJD(3) * t140;
t183 = qJD(1) * t147;
t172 = t142 * t183;
t179 = qJD(3) * t148;
t104 = ((t140 * t179 + t172) * t134 + t179 * t190) * t130;
t193 = t128 * t142;
t99 = (-t104 * t148 + qJD(3)) * t193 + (t172 + (-t104 + t179) * t140) * t127;
t205 = t110 * t111 * t99;
t206 = 0.1e1 / t109 ^ 2 * (-t189 * t205 + (-t142 * t145 * t181 + t164) * t111);
t167 = qJD(1) * t140 + qJD(6);
t180 = qJD(3) * t147;
t156 = t142 * t180 + t167 * t148;
t168 = qJD(6) * t140 + qJD(1);
t161 = t141 * t168;
t105 = t156 * t139 + t147 * t161;
t185 = t148 * t141;
t187 = t147 * t139;
t123 = t140 * t187 - t185;
t119 = t123 ^ 2;
t118 = t119 * t121 + 0.1e1;
t195 = t121 * t123;
t162 = t139 * t168;
t106 = t156 * t141 - t147 * t162;
t201 = t106 * t120 * t121;
t204 = (t105 * t195 - t119 * t201) / t118 ^ 2;
t203 = t104 * t127;
t202 = t104 * t142;
t200 = t111 * t142;
t199 = t111 * t147;
t191 = t134 * t142;
t159 = qJD(3) * (-t134 * t208 - t191);
t198 = (-t135 * t164 + t146 * t159) / t132 ^ 2;
t171 = 0.1e1 + t190;
t117 = t171 * t148 * t130;
t197 = t117 * t148;
t196 = t120 * t139;
t194 = t123 * t141;
t192 = t134 * t138;
t188 = t140 * t148;
t178 = -0.2e1 * t205;
t177 = 0.2e1 * t204;
t176 = t142 * t206;
t175 = t142 * t199;
t174 = t142 * t198;
t173 = t130 * t192;
t170 = 0.2e1 * t123 * t201;
t169 = t198 * t209;
t166 = t148 * t173;
t165 = t207 * t142 * t127;
t163 = t171 * t147;
t160 = t121 * t194 - t196;
t158 = t160 * t147;
t157 = t142 * t179 - t167 * t147;
t126 = t140 * t185 - t187;
t125 = t139 * t188 + t186;
t115 = 0.1e1 / t118;
t107 = 0.1e1 / t109;
t103 = (-t128 * t166 - t165) * t147;
t102 = t127 * t188 + t193 + (-t127 * t140 - t128 * t184) * t117;
t100 = -t171 * t169 + (-qJD(1) * t163 + t159 * t209) * t130;
t1 = [-0.2e1 * t147 * t134 * t174 + (-qJD(3) * t163 + t182 * t191) * t130, 0, t100, 0, 0, 0; (0.2e1 * t110 * t176 + (t110 * t181 + (qJD(1) * t103 + t99) * t200) * t107) * t148 + (-0.2e1 * t111 * t176 * t103 + (((t104 * t166 + t207 * t181 + 0.2e1 * t174) * t127 + (t169 * t192 + t202 + (-t202 + (0.2e1 * t142 + t208) * t179) * t130) * t128) * t175 + (-t111 * t181 + t142 * t178) * t103 + (t110 + ((t145 - t146) * t128 * t173 - t148 * t165) * t111) * t142 * qJD(1)) * t107) * t147, 0, 0.2e1 * (-t102 * t200 - t110 * t140) * t147 * t206 + ((t110 * t182 + (-qJD(3) * t102 - t99) * t199) * t140 + (t110 * t180 + (-t100 * t128 * t148 + t127 * t179 + t197 * t203 - t203 + (-qJD(3) * t127 + t128 * t183) * t117) * t175 + (t111 * t182 + t147 * t178) * t102 + ((-t100 - t183) * t127 + ((-0.1e1 + t197) * qJD(3) + (-t117 + t148) * t104) * t128) * t140 * t199) * t142) * t107, 0, 0, 0; (-t120 * t125 + t126 * t195) * t177 + (t126 * t170 + t148 * t120 * t161 + t157 * t196 + (t148 * t123 * t162 - t126 * t105 - t125 * t106 - t157 * t194) * t121) * t115, 0, t142 * t158 * t177 + (t158 * t181 + (-t160 * t182 + ((qJD(6) * t120 + t170) * t141 + (-t105 * t141 + (qJD(6) * t123 - t106) * t139) * t121) * t147) * t142) * t115, 0, 0, -0.2e1 * t204 + 0.2e1 * (t105 * t115 * t121 + (-t115 * t201 - t121 * t204) * t123) * t123;];
JaD_rot  = t1;
