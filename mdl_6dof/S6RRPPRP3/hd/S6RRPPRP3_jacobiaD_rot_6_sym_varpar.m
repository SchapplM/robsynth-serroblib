% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRP3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_jacobiaD_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:09
% EndTime: 2019-02-26 21:26:10
% DurationCPUTime: 0.73s
% Computational Cost: add. (813->93), mult. (2519->208), div. (480->12), fcn. (2968->9), ass. (0->92)
t134 = sin(qJ(2));
t135 = sin(qJ(1));
t137 = cos(qJ(2));
t181 = t135 * t137;
t119 = atan2(-t181, t134);
t118 = cos(t119);
t117 = sin(t119);
t167 = t117 * t181;
t103 = t118 * t134 - t167;
t100 = 0.1e1 / t103;
t136 = cos(qJ(5));
t138 = cos(qJ(1));
t180 = t136 * t138;
t165 = t134 * t180;
t133 = sin(qJ(5));
t183 = t135 * t133;
t116 = t165 - t183;
t110 = 0.1e1 / t116;
t126 = 0.1e1 / t134;
t101 = 0.1e1 / t103 ^ 2;
t111 = 0.1e1 / t116 ^ 2;
t127 = 0.1e1 / t134 ^ 2;
t129 = t135 ^ 2;
t131 = t137 ^ 2;
t187 = t127 * t131;
t122 = t129 * t187 + 0.1e1;
t120 = 0.1e1 / t122;
t198 = t120 - 0.1e1;
t132 = t138 ^ 2;
t186 = t131 * t132;
t99 = t101 * t186 + 0.1e1;
t97 = 0.1e1 / t99;
t197 = t101 * t97;
t174 = qJD(2) * t135;
t163 = t127 * t174;
t176 = qJD(1) * t138;
t164 = t137 * t176;
t94 = ((t134 * t174 - t164) * t126 + t131 * t163) * t120;
t155 = -t94 + t174;
t156 = -t135 * t94 + qJD(2);
t189 = t118 * t137;
t89 = t156 * t189 + (t155 * t134 - t164) * t117;
t196 = t100 * t101 * t89;
t182 = t135 * t136;
t184 = t134 * t138;
t115 = t133 * t184 + t182;
t109 = t115 ^ 2;
t108 = t109 * t111 + 0.1e1;
t192 = t111 * t115;
t153 = qJD(1) * t134 + qJD(5);
t154 = qJD(5) * t134 + qJD(1);
t173 = qJD(2) * t137;
t162 = t138 * t173;
t185 = t133 * t138;
t96 = -t154 * t185 + (-t153 * t135 + t162) * t136;
t194 = t110 * t111 * t96;
t95 = -qJD(5) * t165 - t133 * t162 - t136 * t176 + t153 * t183;
t195 = 0.1e1 / t108 ^ 2 * (-t109 * t194 - t95 * t192);
t193 = t110 * t133;
t191 = t115 * t136;
t190 = t117 * t134;
t188 = t126 * t131;
t178 = qJD(1) * t135;
t177 = qJD(1) * t137;
t175 = qJD(2) * t134;
t151 = t131 * t135 * t176;
t172 = 0.2e1 * (-t186 * t196 + (-t132 * t134 * t173 - t151) * t101) / t99 ^ 2;
t171 = 0.2e1 * t196;
t170 = 0.2e1 * t195;
t130 = t137 * t131;
t147 = qJD(2) * (-t127 * t130 - t137) * t126;
t169 = 0.2e1 * (t127 * t151 + t129 * t147) / t122 ^ 2;
t168 = t97 * t175;
t166 = t120 * t188;
t161 = 0.1e1 + t187;
t160 = t100 * t172;
t159 = 0.2e1 * t115 * t194;
t158 = t135 * t169;
t157 = t137 * t169;
t152 = t135 * t166;
t150 = t161 * t138;
t149 = t153 * t138;
t148 = t111 * t191 - t193;
t146 = t148 * t138;
t114 = -t134 * t182 - t185;
t113 = -t134 * t183 + t180;
t107 = t161 * t135 * t120;
t105 = 0.1e1 / t108;
t93 = (t198 * t137 * t117 + t118 * t152) * t138;
t92 = t135 * t190 + t189 + (-t118 * t181 - t190) * t107;
t90 = -t161 * t158 + (qJD(1) * t150 + 0.2e1 * t135 * t147) * t120;
t1 = [t126 * t138 * t157 + (t126 * t135 * t177 + qJD(2) * t150) * t120, t90, 0, 0, 0, 0; (t100 * t168 + (t160 + (qJD(1) * t93 + t89) * t197) * t137) * t135 + (t93 * t137 * t97 * t171 + (t93 * t168 + (t93 * t172 + ((t94 * t152 + t198 * t175 + t157) * t117 + (t158 * t188 + t137 * t94 + (t130 * t163 - (t94 - 0.2e1 * t174) * t137) * t120) * t118) * t97 * t138) * t137) * t101 + (-t100 + ((t129 - t132) * t118 * t166 + t198 * t167) * t101) * t97 * t177) * t138 (t100 * t97 * t178 + (t160 + (qJD(2) * t92 + t89) * t197) * t138) * t134 + (t92 * t138 * t101 * t172 + ((-qJD(2) * t100 + t92 * t171) * t138 + (t92 * t178 + (-(-t107 * t176 - t135 * t90) * t118 - ((t107 * t135 - 0.1e1) * t94 + (-t107 + t135) * qJD(2)) * t117) * t137 * t138) * t101) * t97 - ((-t90 + t176) * t117 + (t155 * t107 - t156) * t118) * t184 * t197) * t137, 0, 0, 0, 0; (-t110 * t113 + t114 * t192) * t170 + (t114 * t159 - t154 * t110 * t182 + (-t135 * t173 - t149) * t193 + (-t113 * t96 + t114 * t95 + t149 * t191 - (t154 * t133 - t136 * t173) * t115 * t135) * t111) * t105, t137 * t146 * t170 + (t146 * t175 + (t148 * t178 + ((qJD(5) * t110 + t159) * t136 + (t136 * t95 + (qJD(5) * t115 - t96) * t133) * t111) * t138) * t137) * t105, 0, 0, -0.2e1 * t195 + 0.2e1 * (-t105 * t111 * t95 + (-t105 * t194 - t111 * t195) * t115) * t115, 0;];
JaD_rot  = t1;
