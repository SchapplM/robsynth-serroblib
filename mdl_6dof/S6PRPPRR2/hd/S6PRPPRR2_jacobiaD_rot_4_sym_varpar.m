% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPPRR2_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:22
% EndTime: 2019-02-26 19:45:23
% DurationCPUTime: 0.32s
% Computational Cost: add. (1331->44), mult. (4118->102), div. (242->14), fcn. (5370->11), ass. (0->57)
t122 = cos(pkin(6));
t117 = sin(pkin(11));
t120 = cos(pkin(11));
t123 = sin(qJ(2));
t124 = cos(qJ(2));
t136 = t124 * t117 + t123 * t120;
t110 = t136 * t122;
t113 = t123 * t117 - t124 * t120;
t111 = t113 * qJD(2);
t119 = sin(pkin(6));
t108 = t113 * t119;
t118 = sin(pkin(10));
t121 = cos(pkin(10));
t134 = t113 * t122;
t95 = -t118 * t136 - t121 * t134;
t87 = atan2(t95, t108);
t82 = sin(t87);
t83 = cos(t87);
t81 = t83 * t108 + t82 * t95;
t78 = 0.1e1 / t81;
t97 = t118 * t134 - t121 * t136;
t149 = -0.2e1 * t97;
t105 = 0.1e1 / t108;
t106 = 0.1e1 / t108 ^ 2;
t79 = 0.1e1 / t81 ^ 2;
t109 = t136 * t119;
t102 = qJD(2) * t109;
t139 = -t108 * t82 + t83 * t95;
t146 = t106 * t95;
t92 = t95 ^ 2;
t86 = t92 * t106 + 0.1e1;
t84 = 0.1e1 / t86;
t133 = qJD(2) * t110;
t89 = t118 * t111 - t121 * t133;
t72 = (-t102 * t146 + t105 * t89) * t84;
t70 = t83 * t102 + t139 * t72 + t82 * t89;
t148 = t70 * t78 * t79;
t147 = t79 * t97;
t145 = t109 * t95;
t144 = t102 * t105 * t106;
t140 = 0.1e1 / t118 ^ 2 / t119 ^ 2;
t104 = t122 * t111;
t112 = t136 * qJD(2);
t138 = t118 * t104 - t121 * t112;
t137 = -t118 * t110 - t121 * t113;
t94 = -t121 * t110 + t118 * t113;
t135 = -t105 * t94 + t106 * t145;
t103 = t119 * t111;
t93 = t97 ^ 2;
t91 = t121 * t111 + t118 * t133;
t90 = t121 * t104 + t118 * t112;
t88 = t137 ^ 2 * t140 + 0.1e1;
t76 = t93 * t79 + 0.1e1;
t73 = t135 * t84;
t71 = t83 * t109 - t139 * t73 + t82 * t94;
t69 = 0.2e1 * t135 / t86 ^ 2 * (-t92 * t144 + t89 * t146) + (0.2e1 * t144 * t145 + t105 * t90 + (-t102 * t94 + t103 * t95 - t109 * t89) * t106) * t84;
t1 = [0, t69, 0, 0, 0, 0; 0, 0.2e1 * (-t137 * t78 - t71 * t147) / t76 ^ 2 * (t91 * t147 - t93 * t148) + (t138 * t78 + t71 * t148 * t149 + (-t137 * t70 + t71 * t91 + ((t69 * t95 - t73 * t89 - t103 + (t108 * t73 + t94) * t72) * t83 + (t102 * t73 - t108 * t69 + t90 + (t73 * t95 - t109) * t72) * t82) * t97) * t79) / t76, 0, 0, 0, 0; 0 (t91 / t88 + 0.1e1 / t88 ^ 2 * t137 * t138 * t140 * t149) / t118 / t119, 0, 0, 0, 0;];
JaD_rot  = t1;
