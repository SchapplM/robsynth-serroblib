% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP3_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:09:08
% EndTime: 2019-02-26 21:09:08
% DurationCPUTime: 0.70s
% Computational Cost: add. (1976->91), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->94)
t123 = qJ(1) + pkin(10);
t121 = sin(t123);
t119 = t121 ^ 2;
t130 = sin(qJ(3));
t125 = t130 ^ 2;
t132 = cos(qJ(3));
t127 = 0.1e1 / t132 ^ 2;
t175 = t125 * t127;
t115 = t119 * t175 + 0.1e1;
t124 = t130 * t125;
t126 = 0.1e1 / t132;
t140 = qJD(3) * (t124 * t127 + t130) * t126;
t122 = cos(t123);
t171 = qJD(1) * t122;
t180 = t121 * t125;
t145 = t171 * t180;
t187 = 0.1e1 / t115 ^ 2 * (t119 * t140 + t127 * t145);
t198 = -0.2e1 * t187;
t129 = sin(qJ(4));
t131 = cos(qJ(4));
t173 = t131 * t132;
t109 = t121 * t129 + t122 * t173;
t103 = 0.1e1 / t109;
t104 = 0.1e1 / t109 ^ 2;
t174 = t129 * t132;
t108 = -t121 * t131 + t122 * t174;
t183 = t104 * t108;
t142 = -t103 * t129 + t131 * t183;
t102 = t108 ^ 2;
t101 = t102 * t104 + 0.1e1;
t98 = 0.1e1 / t101;
t197 = t142 * t98;
t169 = qJD(3) * t121;
t113 = 0.1e1 / t115;
t156 = t127 * t169;
t170 = qJD(1) * t130;
t157 = t122 * t170;
t167 = qJD(3) * t132;
t87 = (-(-t121 * t167 - t157) * t126 + t125 * t156) * t113;
t149 = t87 - t169;
t152 = 0.1e1 + t175;
t196 = t121 * t152;
t150 = -t121 * t87 + qJD(3);
t148 = qJD(4) * t132 - qJD(1);
t168 = qJD(3) * t130;
t195 = t129 * t148 + t131 * t168;
t178 = t121 * t130;
t112 = atan2(-t178, -t132);
t111 = cos(t112);
t110 = sin(t112);
t161 = t110 * t178;
t96 = -t111 * t132 - t161;
t93 = 0.1e1 / t96;
t94 = 0.1e1 / t96 ^ 2;
t194 = t113 - 0.1e1;
t120 = t122 ^ 2;
t162 = t94 * t167;
t181 = t111 * t130;
t82 = t150 * t181 + (t132 * t149 - t157) * t110;
t192 = t82 * t93 * t94;
t92 = t120 * t125 * t94 + 0.1e1;
t193 = (-t94 * t145 + (-t125 * t192 + t130 * t162) * t120) / t92 ^ 2;
t147 = -qJD(1) * t132 + qJD(4);
t143 = t147 * t131;
t89 = t121 * t143 - t195 * t122;
t188 = t103 * t104 * t89;
t141 = t121 * t174 + t122 * t131;
t155 = t129 * t168;
t88 = t141 * qJD(1) - t109 * qJD(4) + t122 * t155;
t191 = (-t102 * t188 - t183 * t88) / t101 ^ 2;
t90 = 0.1e1 / t92;
t190 = t90 * t94;
t189 = t93 * t90;
t185 = t122 * t94;
t182 = t110 * t132;
t177 = t122 * t129;
t176 = t122 * t130;
t172 = qJD(1) * t121;
t166 = 0.2e1 * t192;
t165 = -0.2e1 * t191;
t164 = t93 * t193;
t163 = t108 * t188;
t160 = t113 * t125 * t126;
t158 = t121 * t170;
t153 = 0.2e1 * t94 * t193;
t151 = t126 * t198;
t146 = t121 * t160;
t144 = t152 * t122;
t107 = -t121 * t173 + t177;
t100 = t113 * t196;
t86 = (t110 * t130 * t194 - t111 * t146) * t122;
t85 = -t121 * t182 + t181 + (-t111 * t178 + t182) * t100;
t83 = t196 * t198 + (qJD(1) * t144 + 0.2e1 * t121 * t140) * t113;
t1 = [t151 * t176 + (qJD(3) * t144 - t126 * t158) * t113, 0, t83, 0, 0, 0; (-t167 * t189 + (0.2e1 * t164 + (qJD(1) * t86 + t82) * t190) * t130) * t121 + (t86 * t153 * t130 + (-t86 * t162 + (t86 * t166 + ((0.2e1 * t130 * t187 - t146 * t87 - t167 * t194) * t110 + (t151 * t180 + t130 * t87 + (t124 * t156 - (t87 - 0.2e1 * t169) * t130) * t113) * t111) * t185) * t130 + (-t93 + (-(t119 - t120) * t111 * t160 + t194 * t161) * t94) * t170) * t90) * t122, 0 (-t172 * t189 + (-0.2e1 * t164 + (-qJD(3) * t85 - t82) * t190) * t122) * t132 + (t85 * t122 * t153 + (-t122 * qJD(3) * t93 - ((-t100 * t171 - t121 * t83) * t111 + (-t150 * t100 - t149) * t110) * t94 * t176 + (t122 * t166 + t172 * t94) * t85 - ((t83 - t171) * t110 + (t100 * t149 + t150) * t111) * t132 * t185) * t90) * t130, 0, 0, 0; 0.2e1 * (t103 * t141 + t107 * t183) * t191 + (0.2e1 * t107 * t163 + (t141 * t89 + t107 * t88 + (-t195 * t121 - t122 * t143) * t108) * t104 + (t147 * t177 + (-t131 * t148 + t155) * t121) * t103) * t98, 0, -t158 * t197 + (t167 * t197 + (t142 * t165 + ((-qJD(4) * t103 - 0.2e1 * t163) * t131 + (-t131 * t88 + (-qJD(4) * t108 + t89) * t129) * t104) * t98) * t130) * t122, t165 + 0.2e1 * (-t104 * t88 * t98 + (-t104 * t191 - t188 * t98) * t108) * t108, 0, 0;];
JaD_rot  = t1;
