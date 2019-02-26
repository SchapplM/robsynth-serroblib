% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRRPRR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR5_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_jacobiaD_rot_3_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:06:27
% EndTime: 2019-02-26 20:06:28
% DurationCPUTime: 0.38s
% Computational Cost: add. (756->55), mult. (2271->133), div. (423->14), fcn. (2956->11), ass. (0->65)
t138 = sin(qJ(2));
t140 = cos(qJ(2));
t135 = sin(pkin(11));
t165 = cos(pkin(6));
t154 = t135 * t165;
t164 = cos(pkin(11));
t126 = -t138 * t154 + t164 * t140;
t137 = sin(qJ(3));
t139 = cos(qJ(3));
t136 = sin(pkin(6));
t159 = t135 * t136;
t148 = -t126 * t137 + t139 * t159;
t169 = t148 * qJD(3);
t151 = t165 * t164;
t122 = t135 * t138 - t140 * t151;
t158 = t136 * t140;
t112 = atan2(-t122, -t158);
t110 = sin(t112);
t111 = cos(t112);
t97 = -t110 * t122 - t111 * t158;
t94 = 0.1e1 / t97;
t109 = t126 * t139 + t137 * t159;
t105 = 0.1e1 / t109;
t132 = 0.1e1 / t140;
t106 = 0.1e1 / t109 ^ 2;
t133 = 0.1e1 / t140 ^ 2;
t95 = 0.1e1 / t97 ^ 2;
t104 = t148 ^ 2;
t101 = t104 * t106 + 0.1e1;
t147 = -t164 * t138 - t140 * t154;
t118 = t147 * qJD(2);
t102 = t109 * qJD(3) + t118 * t137;
t162 = t106 * t148;
t103 = t118 * t139 + t169;
t163 = t103 * t105 * t106;
t168 = 0.1e1 / t101 ^ 2 * (-t102 * t162 - t104 * t163);
t124 = t135 * t140 + t138 * t151;
t160 = t133 * t138;
t155 = t122 * t160;
t149 = t124 * t132 + t155;
t120 = t122 ^ 2;
t131 = 0.1e1 / t136 ^ 2;
t115 = t120 * t131 * t133 + 0.1e1;
t113 = 0.1e1 / t115;
t130 = 0.1e1 / t136;
t161 = t113 * t130;
t90 = t149 * t161;
t167 = t122 * t90;
t166 = t147 * t95;
t157 = qJD(2) * t138;
t156 = -0.2e1 * t168;
t150 = -t105 * t137 - t139 * t162;
t134 = t132 * t133;
t121 = t147 ^ 2;
t119 = t126 * qJD(2);
t117 = t124 * qJD(2);
t116 = t122 * qJD(2);
t99 = 0.1e1 / t101;
t96 = t94 * t95;
t93 = t121 * t95 + 0.1e1;
t89 = (qJD(2) * t155 + t117 * t132) * t161;
t87 = (t136 * t138 - t167) * t111 + (t90 * t158 - t124) * t110;
t86 = (-t122 * t89 + t136 * t157) * t111 + (t89 * t158 - t117) * t110;
t85 = (-0.2e1 * t149 * (t117 * t122 * t133 + t120 * t134 * t157) * t131 / t115 ^ 2 + (t117 * t160 - t116 * t132 + (t124 * t160 + (0.2e1 * t134 * t138 ^ 2 + t132) * t122) * qJD(2)) * t113) * t130;
t1 = [0, t85, 0, 0, 0, 0; 0, 0.2e1 * (-t126 * t94 - t87 * t166) / t93 ^ 2 * (-t121 * t96 * t86 - t119 * t166) + (-t87 * t119 * t95 + t118 * t94 + (-0.2e1 * t147 * t87 * t96 - t126 * t95) * t86 - (-(-t117 * t90 - t122 * t85 - t124 * t89 + (t89 * t90 + qJD(2)) * t158) * t111 - (t89 * t167 + t116 + (t140 * t85 + (-qJD(2) * t90 - t89) * t138) * t136) * t110) * t166) / t93, 0, 0, 0, 0; 0, t150 * t99 * t119 - (t150 * t156 + ((-qJD(3) * t105 + 0.2e1 * t148 * t163) * t139 + (t102 * t139 + (t103 + t169) * t137) * t106) * t99) * t147, t156 - 0.2e1 * (t102 * t106 * t99 - (-t106 * t168 - t99 * t163) * t148) * t148, 0, 0, 0;];
JaD_rot  = t1;
