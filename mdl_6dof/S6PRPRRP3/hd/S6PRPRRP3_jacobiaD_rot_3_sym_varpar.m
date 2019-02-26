% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRP3_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_jacobiaD_rot_3_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:33
% EndTime: 2019-02-26 19:51:34
% DurationCPUTime: 0.34s
% Computational Cost: add. (610->48), mult. (1820->126), div. (396->14), fcn. (2421->11), ass. (0->64)
t116 = sin(qJ(2));
t117 = cos(qJ(2));
t113 = sin(pkin(10));
t136 = cos(pkin(6));
t129 = t113 * t136;
t135 = cos(pkin(10));
t104 = -t116 * t129 + t135 * t117;
t126 = t136 * t135;
t100 = t113 * t116 - t117 * t126;
t114 = sin(pkin(6));
t132 = t114 * t117;
t90 = atan2(-t100, -t132);
t88 = sin(t90);
t89 = cos(t90);
t77 = -t88 * t100 - t89 * t132;
t74 = 0.1e1 / t77;
t112 = sin(pkin(11));
t115 = cos(pkin(11));
t133 = t113 * t114;
t87 = t104 * t115 + t112 * t133;
t83 = 0.1e1 / t87;
t109 = 0.1e1 / t117;
t110 = 0.1e1 / t117 ^ 2;
t75 = 0.1e1 / t77 ^ 2;
t84 = 0.1e1 / t87 ^ 2;
t131 = qJD(2) * t116;
t137 = t117 * t88;
t143 = t100 * t89;
t134 = t110 * t116;
t130 = t100 * t134;
t107 = 0.1e1 / t114;
t108 = 0.1e1 / t114 ^ 2;
t98 = t100 ^ 2;
t93 = t98 * t108 * t110 + 0.1e1;
t91 = 0.1e1 / t93;
t141 = t107 * t91;
t102 = t113 * t117 + t116 * t126;
t95 = t102 * qJD(2);
t69 = (qJD(2) * t130 + t109 * t95) * t141;
t66 = -t69 * t143 - t88 * t95 + (t89 * t131 + t69 * t137) * t114;
t146 = t66 * t74 * t75;
t139 = t112 * t84;
t86 = t104 * t112 - t115 * t133;
t82 = t86 ^ 2;
t81 = t82 * t84 + 0.1e1;
t85 = t83 * t84;
t124 = -t135 * t116 - t117 * t129;
t96 = t124 * qJD(2);
t145 = (-t115 * t82 * t85 + t86 * t139) * t96 / t81 ^ 2;
t125 = t102 * t109 + t130;
t70 = t125 * t141;
t144 = t69 * t70;
t142 = t124 * t75;
t140 = t112 * t83;
t138 = t115 * t86;
t111 = t109 * t110;
t99 = t124 ^ 2;
t97 = t104 * qJD(2);
t94 = t100 * qJD(2);
t79 = 0.1e1 / t81;
t73 = t99 * t75 + 0.1e1;
t67 = -t70 * t143 - t88 * t102 + (t116 * t89 + t70 * t137) * t114;
t65 = (-0.2e1 * t125 / t93 ^ 2 * (t100 * t110 * t95 + t111 * t98 * t131) * t108 + (t95 * t134 - t109 * t94 + (t102 * t134 + (0.2e1 * t111 * t116 ^ 2 + t109) * t100) * qJD(2)) * t91) * t107;
t1 = [0, t65, 0, 0, 0, 0; 0, 0.2e1 * (-t104 * t74 - t67 * t142) / t73 ^ 2 * (-t97 * t142 - t99 * t146) + (-0.2e1 * t67 * t124 * t146 + t96 * t74 + (-t104 * t66 - t67 * t97 - (-(t100 * t144 + t94) * t88 - (-t100 * t65 - t102 * t69 - t70 * t95) * t89) * t124) * t75 + ((-qJD(2) * t70 - t69) * t88 * t116 + (t65 * t88 + (qJD(2) + t144) * t89) * t117) * t114 * t142) / t73, 0, 0, 0, 0; 0 (t84 * t138 - t140) * t97 * t79 - 0.2e1 * (t140 * t145 + (-t84 * t86 * t145 + (-t85 * t138 + t139) * t96 * t79) * t115) * t124, 0, 0, 0, 0;];
JaD_rot  = t1;
