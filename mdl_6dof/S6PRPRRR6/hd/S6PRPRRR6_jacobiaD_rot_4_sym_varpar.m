% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPRRR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR6_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:50
% EndTime: 2019-02-26 19:56:51
% DurationCPUTime: 0.36s
% Computational Cost: add. (682->54), mult. (2271->134), div. (423->14), fcn. (2956->11), ass. (0->64)
t135 = sin(pkin(11));
t137 = cos(pkin(11));
t140 = sin(qJ(2));
t138 = cos(pkin(6));
t142 = cos(qJ(2));
t155 = t138 * t142;
t123 = t135 * t140 - t137 * t155;
t156 = t138 * t140;
t124 = t135 * t142 + t137 * t156;
t136 = sin(pkin(6));
t157 = t136 * t140;
t114 = atan2(-t124, t157);
t110 = sin(t114);
t111 = cos(t114);
t97 = -t110 * t124 + t111 * t157;
t94 = 0.1e1 / t97;
t126 = t135 * t155 + t137 * t140;
t139 = sin(qJ(4));
t141 = cos(qJ(4));
t159 = t135 * t136;
t109 = t126 * t139 + t141 * t159;
t105 = 0.1e1 / t109;
t132 = 0.1e1 / t140;
t106 = 0.1e1 / t109 ^ 2;
t133 = 0.1e1 / t140 ^ 2;
t95 = 0.1e1 / t97 ^ 2;
t108 = -t126 * t141 + t139 * t159;
t104 = t108 ^ 2;
t101 = t104 * t106 + 0.1e1;
t127 = -t135 * t156 + t137 * t142;
t120 = t127 * qJD(2);
t103 = t109 * qJD(4) - t120 * t141;
t162 = t106 * t108;
t153 = qJD(4) * t108;
t102 = t120 * t139 - t153;
t163 = t102 * t105 * t106;
t166 = 0.1e1 / t101 ^ 2 * (t103 * t162 - t104 * t163);
t160 = t133 * t142;
t152 = t124 * t160;
t149 = t123 * t132 + t152;
t121 = t124 ^ 2;
t131 = 0.1e1 / t136 ^ 2;
t115 = t121 * t131 * t133 + 0.1e1;
t112 = 0.1e1 / t115;
t130 = 0.1e1 / t136;
t161 = t112 * t130;
t90 = t149 * t161;
t165 = t124 * t90;
t164 = t127 * t95;
t154 = qJD(2) * t142;
t150 = t105 * t141 + t139 * t162;
t134 = t132 * t133;
t122 = t127 ^ 2;
t119 = t126 * qJD(2);
t118 = t124 * qJD(2);
t117 = t123 * qJD(2);
t99 = 0.1e1 / t101;
t96 = t94 * t95;
t93 = t122 * t95 + 0.1e1;
t89 = (qJD(2) * t152 + t117 * t132) * t161;
t87 = (t136 * t142 - t165) * t111 + (-t90 * t157 + t123) * t110;
t86 = (-t124 * t89 + t136 * t154) * t111 + (-t89 * t157 + t117) * t110;
t85 = (-0.2e1 * t149 * (-t117 * t124 * t133 - t121 * t134 * t154) * t131 / t115 ^ 2 + (-t117 * t160 + t118 * t132 + (-t123 * t160 + (-0.2e1 * t134 * t142 ^ 2 - t132) * t124) * qJD(2)) * t112) * t130;
t1 = [0, t85, 0, 0, 0, 0; 0, 0.2e1 * (t126 * t94 + t87 * t164) / t93 ^ 2 * (-t122 * t86 * t96 - t119 * t164) + (t87 * t119 * t95 - t120 * t94 + (0.2e1 * t127 * t87 * t96 + t126 * t95) * t86 + (-(t117 * t90 + t123 * t89 - t124 * t85 + (-t89 * t90 - qJD(2)) * t157) * t111 - (t89 * t165 + t118 + (-t140 * t85 + (-qJD(2) * t90 - t89) * t142) * t136) * t110) * t164) / t93, 0, 0, 0, 0; 0, t150 * t99 * t119 + (0.2e1 * t150 * t166 + ((qJD(4) * t105 + 0.2e1 * t108 * t163) * t139 + (-t103 * t139 + (t102 - t153) * t141) * t106) * t99) * t127, 0, -0.2e1 * t166 + 0.2e1 * (t103 * t106 * t99 + (-t106 * t166 - t99 * t163) * t108) * t108, 0, 0;];
JaD_rot  = t1;
