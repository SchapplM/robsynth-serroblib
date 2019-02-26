% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
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

function JaD_rot = S4RPPP1_jacobiaD_rot_3_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobiaD_rot_3_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_jacobiaD_rot_3_floatb_twist_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobiaD_rot_3_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:47
% EndTime: 2018-11-14 13:45:47
% DurationCPUTime: 0.24s
% Computational Cost: add. (1267->36), mult. (1366->90), div. (132->14), fcn. (1422->12), ass. (0->51)
t105 = cos(qJ(1));
t104 = sin(qJ(1));
t99 = 0.1e1 / t104 ^ 2;
t120 = t105 * t99;
t101 = sin(pkin(6));
t116 = pkin(4) + pkin(6);
t117 = pkin(4) - pkin(6);
t113 = cos(t117) / 0.2e1 + cos(t116) / 0.2e1;
t82 = t104 * t101 - t105 * t113;
t94 = -sin(t117) / 0.2e1;
t95 = sin(t116);
t93 = -t95 / 0.2e1 + t94;
t73 = atan2(-t82, t93);
t68 = sin(t73);
t69 = cos(t73);
t127 = t82 ^ 2;
t90 = 0.1e1 / t93 ^ 2;
t72 = t127 * t90 + 0.1e1;
t70 = 0.1e1 / t72;
t89 = 0.1e1 / t93;
t114 = -t68 + (t69 * t82 * t89 + t68) * t70;
t66 = -t68 * t82 + t69 * t93;
t63 = 0.1e1 / t66;
t98 = 0.1e1 / t104;
t64 = 0.1e1 / t66 ^ 2;
t85 = t105 * t101 + t104 * t113;
t79 = t85 * qJD(1);
t58 = t114 * t79;
t126 = t58 * t63 * t64;
t125 = t70 * t89;
t71 = 0.1e1 / t72 ^ 2;
t124 = t71 * t82;
t77 = t82 * qJD(1);
t123 = t77 * t64;
t122 = t79 * t85;
t92 = t95 / 0.2e1 + t94;
t121 = t104 * t92;
t119 = t98 * t120;
t103 = cos(pkin(6));
t84 = -t104 * t103 - t105 * t92;
t102 = sin(pkin(4));
t97 = 0.1e1 / t102 ^ 2;
t91 = t89 * t90;
t86 = t105 * t103 - t121;
t81 = t86 ^ 2;
t80 = t85 ^ 2;
t78 = t84 * qJD(1);
t76 = t81 * t99 * t97 + 0.1e1;
t62 = t80 * t64 + 0.1e1;
t59 = t114 * t85;
t1 = [0.2e1 * t91 * t122 * t124 + t77 * t125, 0, 0, 0; 0.2e1 * (t59 * t64 * t85 + t63 * t82) / t62 ^ 2 * (-t85 * t123 - t80 * t126) + (-t79 * t63 + (t82 * t58 + t59 * t77) * t64 + (0.2e1 * t59 * t126 + t114 * t123 - (-t68 * t90 * t124 + (0.2e1 * t125 + (-0.2e1 * t127 * t91 - t89) * t71) * t69) * t64 * t122) * t85) / t62, 0, 0, 0; (0.2e1 * (t86 * t120 - t84 * t98) / t76 ^ 2 * (-qJD(1) * t81 * t119 + t78 * t86 * t99) * t97 + (-t78 * t120 + ((t86 + t121) * t98 + (-t103 * t98 + 0.2e1 * t86 * t119 - t84 * t99) * t105) * qJD(1)) / t76) / t102, 0, 0, 0;];
JaD_rot  = t1;
