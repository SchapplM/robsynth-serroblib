% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR4_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobiaD_rot_3_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:46
% EndTime: 2019-02-26 21:55:46
% DurationCPUTime: 0.36s
% Computational Cost: add. (454->46), mult. (1493->122), div. (133->12), fcn. (1882->11), ass. (0->62)
t133 = sin(qJ(1));
t134 = cos(qJ(1));
t131 = cos(pkin(6));
t129 = sin(pkin(12));
t132 = sin(qJ(2));
t161 = cos(pkin(12));
t166 = cos(qJ(2));
t145 = -t132 * t129 + t166 * t161;
t144 = t145 * t131;
t146 = t166 * t129 + t132 * t161;
t99 = -t133 * t144 - t134 * t146;
t92 = t99 ^ 2;
t109 = t146 * t131;
t147 = t133 * t109 - t134 * t145;
t94 = 0.1e1 / t147 ^ 2;
t169 = t92 * t94;
t125 = 0.1e1 / t131 ^ 2;
t130 = sin(pkin(6));
t123 = t130 ^ 2;
t128 = t134 ^ 2;
t119 = t128 * t123 * t125 + 0.1e1;
t127 = t133 ^ 2;
t159 = 0.1e1 / t119 ^ 2 * t127;
t168 = t125 * t159;
t167 = qJD(1) * t144 + t145 * qJD(2);
t93 = 0.1e1 / t147;
t155 = t134 * t130;
t118 = atan2(t155, t131);
t114 = sin(t118);
t115 = cos(t118);
t105 = t114 * t155 + t115 * t131;
t102 = 0.1e1 / t105;
t124 = 0.1e1 / t131;
t103 = 0.1e1 / t105 ^ 2;
t163 = t94 * t99;
t111 = t146 * qJD(2);
t108 = t131 * t111;
t156 = t133 * t146;
t83 = -qJD(1) * t156 - t133 * t108 + t167 * t134;
t151 = t83 * t163;
t107 = qJD(2) * t144;
t97 = -t134 * t109 - t133 * t145;
t84 = t97 * qJD(1) - t133 * t107 - t134 * t111;
t95 = t93 * t94;
t162 = t95 * t84;
t88 = 0.1e1 + t169;
t165 = (t92 * t162 - t151) / t88 ^ 2;
t160 = t103 * t133;
t158 = t123 * t124;
t154 = qJD(1) * t134;
t152 = -0.2e1 * t97 * t99;
t116 = 0.1e1 / t119;
t150 = (t116 - 0.1e1) * t130;
t149 = -0.2e1 * t124 * t168;
t85 = (-t115 * t116 * t134 * t158 + t114 * t150) * t133;
t122 = t130 * t123;
t104 = t102 * t103;
t96 = t134 * t144 - t156;
t91 = t127 * t123 * t103 + 0.1e1;
t86 = 0.1e1 / t88;
t82 = qJD(1) * t85;
t1 = [(-t116 * t124 * t130 + t122 * t149) * t154, 0, 0, 0, 0, 0; (0.2e1 * (-t102 * t134 + t85 * t160) / t91 ^ 2 * (-t104 * t127 * t82 + t154 * t160) * t123 + ((0.2e1 * t104 * t133 * t85 - t103 * t134) * t82 + (-t133 * t102 + ((-t85 + (-t122 * t168 - t150) * t133 * t114) * t134 - (t128 * t123 ^ 2 * t149 + (-t159 + (0.2e1 * t127 - t128) * t116) * t158) * t133 * t115) * t103) * qJD(1)) / t91) * t130, 0, 0, 0, 0, 0; (t94 * t152 + 0.2e1 * t96 * t93) * t165 + ((-t97 * t83 - t96 * t84) * t94 - (-t134 * t108 - t167 * t133 - t146 * t154) * t93 + (t147 * qJD(1) - t134 * t107 + t133 * t111) * t163 - t152 * t162) * t86, 0.2e1 * (-t147 * t93 - t169) * t165 + (-0.2e1 * t151 + (t147 * t94 + 0.2e1 * t92 * t95 - t93) * t84) * t86, 0, 0, 0, 0;];
JaD_rot  = t1;
