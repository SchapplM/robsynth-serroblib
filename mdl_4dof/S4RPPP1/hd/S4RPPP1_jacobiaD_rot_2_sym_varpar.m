% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4RPPP1_jacobiaD_rot_2_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobiaD_rot_2_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_jacobiaD_rot_2_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobiaD_rot_2_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:29:48
% EndTime: 2019-02-26 19:29:48
% DurationCPUTime: 0.21s
% Computational Cost: add. (137->30), mult. (614->95), div. (108->12), fcn. (792->9), ass. (0->49)
t86 = sin(pkin(4));
t79 = t86 ^ 2;
t88 = cos(pkin(4));
t81 = 0.1e1 / t88 ^ 2;
t90 = cos(qJ(1));
t84 = t90 ^ 2;
t77 = t79 * t81 * t84 + 0.1e1;
t89 = sin(qJ(1));
t83 = t89 ^ 2;
t108 = 0.1e1 / t77 ^ 2 * t83;
t112 = t108 * t81;
t103 = t90 * t86;
t76 = atan2(t103, t88);
t72 = sin(t76);
t73 = cos(t76);
t58 = t72 * t103 + t73 * t88;
t55 = 0.1e1 / t58;
t85 = sin(pkin(6));
t105 = t89 * t85;
t87 = cos(pkin(6));
t99 = t88 * t105 - t87 * t90;
t65 = 0.1e1 / t99;
t80 = 0.1e1 / t88;
t56 = 0.1e1 / t58 ^ 2;
t66 = 0.1e1 / t99 ^ 2;
t111 = t56 * t89;
t104 = t89 * t87;
t70 = t88 * t104 + t85 * t90;
t110 = t66 * t70;
t106 = t88 * t90;
t69 = -t85 * t106 - t104;
t109 = t69 * t70;
t107 = t79 * t80;
t102 = qJD(1) * t90;
t74 = 0.1e1 / t77;
t101 = (t74 - 0.1e1) * t86;
t100 = -0.2e1 * t80 * t112;
t68 = t87 * t106 - t105;
t51 = (-t73 * t74 * t90 * t107 + t72 * t101) * t89;
t78 = t86 * t79;
t67 = t65 * t66;
t64 = t70 ^ 2;
t63 = t69 * qJD(1);
t62 = t68 * qJD(1);
t61 = t64 * t66 + 0.1e1;
t57 = t55 * t56;
t54 = t56 * t79 * t83 + 0.1e1;
t50 = qJD(1) * t51;
t1 = [(-t74 * t80 * t86 + t78 * t100) * t102, 0, 0, 0; (0.2e1 * (t51 * t111 - t55 * t90) / t54 ^ 2 * (-t50 * t57 * t83 + t102 * t111) * t79 + ((0.2e1 * t51 * t57 * t89 - t56 * t90) * t50 + (-t89 * t55 + ((-t51 + (-t78 * t112 - t101) * t89 * t72) * t90 - (t84 * t79 ^ 2 * t100 + (-t108 + (0.2e1 * t83 - t84) * t74) * t107) * t89 * t73) * t56) * qJD(1)) / t54) * t86, 0, 0, 0; 0.2e1 * (t66 * t109 + t65 * t68) / t61 ^ 2 * (t63 * t64 * t67 + t62 * t110) + (-t69 * t62 * t66 + (-0.2e1 * t67 * t109 - t68 * t66) * t63 + (-t99 * t110 + t70 * t65) * qJD(1)) / t61, 0, 0, 0;];
JaD_rot  = t1;
