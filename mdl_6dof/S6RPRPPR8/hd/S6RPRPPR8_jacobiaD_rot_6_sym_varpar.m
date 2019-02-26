% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_jacobiaD_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:43:07
% EndTime: 2019-02-26 20:43:08
% DurationCPUTime: 0.76s
% Computational Cost: add. (624->88), mult. (2519->204), div. (480->12), fcn. (2968->9), ass. (0->91)
t128 = sin(qJ(6));
t131 = cos(qJ(6));
t133 = cos(qJ(1));
t171 = t133 * t131;
t130 = sin(qJ(1));
t132 = cos(qJ(3));
t174 = t130 * t132;
t108 = t128 * t174 - t171;
t196 = 0.2e1 * t108;
t127 = t133 ^ 2;
t129 = sin(qJ(3));
t122 = t129 ^ 2;
t125 = 0.1e1 / t132 ^ 2;
t177 = t122 * t125;
t118 = t127 * t177 + 0.1e1;
t121 = t129 * t122;
t124 = 0.1e1 / t132;
t176 = t124 * t129;
t142 = qJD(3) * (t121 * t124 * t125 + t176);
t169 = qJD(1) * t133;
t160 = t130 * t169;
t185 = 0.1e1 / t118 ^ 2 * (t127 * t142 - t160 * t177);
t195 = -0.2e1 * t185;
t166 = qJD(3) * t133;
t115 = 0.1e1 / t118;
t159 = t125 * t166;
t170 = qJD(1) * t130;
t161 = t129 * t170;
t89 = ((t132 * t166 - t161) * t124 + t122 * t159) * t115;
t152 = -t89 + t166;
t156 = 0.1e1 + t177;
t194 = t133 * t156;
t193 = -t133 * t89 + qJD(3);
t172 = t133 * t129;
t117 = atan2(t172, t132);
t113 = sin(t117);
t114 = cos(t117);
t98 = t113 * t172 + t114 * t132;
t95 = 0.1e1 / t98;
t144 = t128 * t133 + t131 * t174;
t105 = 0.1e1 / t144;
t106 = 0.1e1 / t144 ^ 2;
t96 = 0.1e1 / t98 ^ 2;
t192 = t115 - 0.1e1;
t123 = t130 ^ 2;
t167 = qJD(3) * t132;
t163 = t96 * t167;
t179 = t114 * t129;
t84 = -t193 * t179 + (t152 * t132 - t161) * t113;
t190 = t84 * t95 * t96;
t94 = t122 * t123 * t96 + 0.1e1;
t191 = (t123 * t129 * t163 + (-t123 * t190 + t96 * t160) * t122) / t94 ^ 2;
t92 = 0.1e1 / t94;
t189 = t92 * t96;
t188 = t95 * t92;
t104 = t108 ^ 2;
t103 = t104 * t106 + 0.1e1;
t181 = t106 * t108;
t151 = qJD(6) * t132 + qJD(1);
t146 = t151 * t128;
t150 = qJD(1) * t132 + qJD(6);
t91 = -t150 * t171 + (qJD(3) * t129 * t131 + t146) * t130;
t186 = t105 * t106 * t91;
t173 = t132 * t133;
t143 = t128 * t173 + t130 * t131;
t168 = qJD(3) * t130;
t90 = -t129 * t128 * t168 + t143 * qJD(1) + t144 * qJD(6);
t187 = 0.1e1 / t103 ^ 2 * (t104 * t186 + t90 * t181);
t182 = t105 * t128;
t180 = t108 * t131;
t178 = t122 * t124;
t175 = t129 * t130;
t165 = 0.2e1 * t190;
t164 = -0.2e1 * t187;
t162 = t133 * t178;
t158 = -0.2e1 * t95 * t191;
t157 = 0.2e1 * t96 * t191;
t155 = t186 * t196;
t154 = 0.2e1 * t129 * t185;
t149 = t115 * t162;
t148 = t192 * t129 * t113;
t147 = t156 * t130;
t145 = t106 * t180 - t182;
t141 = t129 * t166 + t150 * t130;
t112 = t130 * t128 - t132 * t171;
t102 = t115 * t194;
t100 = 0.1e1 / t103;
t88 = (-t114 * t149 + t148) * t130;
t87 = t113 * t173 - t179 + (-t113 * t132 + t114 * t172) * t102;
t85 = t194 * t195 + (-qJD(1) * t147 + 0.2e1 * t133 * t142) * t115;
t1 = [t130 * t124 * t154 + (-qJD(3) * t147 - t169 * t176) * t115, 0, t85, 0, 0, 0; (t167 * t188 + (t158 + (-qJD(1) * t88 - t84) * t189) * t129) * t133 + (t88 * t157 * t129 + (-t88 * t163 + (t88 * t165 + ((-t89 * t149 - t192 * t167 + t154) * t113 + (t162 * t195 + t129 * t89 + (t121 * t159 - (t89 - 0.2e1 * t166) * t129) * t115) * t114) * t96 * t130) * t129 + (-t95 + (-(t123 - t127) * t115 * t114 * t178 - t133 * t148) * t96) * t129 * qJD(1)) * t92) * t130, 0 (t169 * t188 + (t158 + (-qJD(3) * t87 - t84) * t189) * t130) * t132 + (t87 * t130 * t157 + (-t95 * t168 - ((-t102 * t170 + t133 * t85) * t114 + (t193 * t102 - t152) * t113) * t96 * t175 + (t130 * t165 - t96 * t169) * t87) * t92 - ((-t85 - t170) * t113 + (t152 * t102 - t193) * t114) * t174 * t189) * t129, 0, 0, 0; 0.2e1 * (-t105 * t143 - t112 * t181) * t187 + (t112 * t155 + t151 * t105 * t171 - t141 * t182 + (t133 * t108 * t146 + t112 * t90 + t141 * t180 + t143 * t91) * t106) * t100, 0, t145 * t164 * t175 + (t145 * t130 * t167 + (t145 * t169 + ((-qJD(6) * t105 + t155) * t131 + (t131 * t90 + (-qJD(6) * t108 - t91) * t128) * t106) * t130) * t129) * t100, 0, 0, t164 + (t100 * t106 * t90 + (t100 * t186 - t106 * t187) * t108) * t196;];
JaD_rot  = t1;
