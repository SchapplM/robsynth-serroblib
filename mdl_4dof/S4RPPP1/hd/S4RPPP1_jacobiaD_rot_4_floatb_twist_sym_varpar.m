% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RPPP1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
%
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S4RPPP1_jacobiaD_rot_4_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobiaD_rot_4_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_jacobiaD_rot_4_floatb_twist_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobiaD_rot_4_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:52
% EndTime: 2018-11-14 13:45:52
% DurationCPUTime: 0.24s
% Computational Cost: add. (1267->34), mult. (1366->90), div. (132->14), fcn. (1422->12), ass. (0->52)
t105 = cos(pkin(6));
t106 = sin(qJ(1));
t107 = cos(qJ(1));
t118 = pkin(4) + pkin(6);
t119 = pkin(4) - pkin(6);
t115 = sin(t118) / 0.2e1 - sin(t119) / 0.2e1;
t85 = t106 * t105 + t107 * t115;
t96 = cos(t119) / 0.2e1;
t97 = cos(t118);
t95 = t96 - t97 / 0.2e1;
t75 = atan2(-t85, t95);
t70 = sin(t75);
t71 = cos(t75);
t127 = t85 ^ 2;
t92 = 0.1e1 / t95 ^ 2;
t74 = t127 * t92 + 0.1e1;
t72 = 0.1e1 / t74;
t91 = 0.1e1 / t95;
t116 = -t70 + (t71 * t85 * t91 + t70) * t72;
t68 = -t70 * t85 + t71 * t95;
t65 = 0.1e1 / t68;
t100 = 0.1e1 / t106;
t101 = 0.1e1 / t106 ^ 2;
t66 = 0.1e1 / t68 ^ 2;
t114 = t106 * t115;
t120 = qJD(1) * t107;
t81 = -qJD(1) * t114 + t105 * t120;
t60 = t116 * t81;
t126 = t60 * t65 * t66;
t125 = t72 * t91;
t73 = 0.1e1 / t74 ^ 2;
t124 = t73 * t85;
t80 = t85 * qJD(1);
t123 = t80 * t66;
t88 = t107 * t105 - t114;
t122 = t81 * t88;
t121 = t101 * t107;
t103 = sin(pkin(6));
t94 = t96 + t97 / 0.2e1;
t87 = -t107 * t103 - t106 * t94;
t84 = t106 * t103 - t107 * t94;
t104 = sin(pkin(4));
t102 = t100 * t101;
t99 = 0.1e1 / t104 ^ 2;
t93 = t91 * t92;
t83 = t88 ^ 2;
t82 = t87 ^ 2;
t79 = t84 * qJD(1);
t78 = t82 * t101 * t99 + 0.1e1;
t64 = t83 * t66 + 0.1e1;
t61 = t116 * t88;
t1 = [0.2e1 * t93 * t122 * t124 + t80 * t125, 0, 0, 0; 0.2e1 * (t61 * t66 * t88 + t65 * t85) / t64 ^ 2 * (-t88 * t123 - t83 * t126) + (-t81 * t65 + (t85 * t60 + t61 * t80) * t66 + (0.2e1 * t61 * t126 + t116 * t123 - (-t70 * t92 * t124 + (0.2e1 * t125 + (-0.2e1 * t127 * t93 - t91) * t73) * t71) * t66 * t122) * t88) / t64, 0, 0, 0; (0.2e1 * (-t100 * t84 + t87 * t121) / t78 ^ 2 * (t101 * t79 * t87 - t102 * t82 * t120) * t99 + (-t79 * t121 + (0.2e1 * t87 * t102 * t107 - t84 * t101) * t120) / t78) / t104, 0, 0, 0;];
JaD_rot  = t1;
