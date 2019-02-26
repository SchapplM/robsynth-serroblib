% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR12_jacobiaD_rot_2_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobiaD_rot_2_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_jacobiaD_rot_2_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobiaD_rot_2_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:07:16
% EndTime: 2019-02-26 21:07:16
% DurationCPUTime: 0.24s
% Computational Cost: add. (137->30), mult. (614->94), div. (108->12), fcn. (792->9), ass. (0->50)
t86 = sin(pkin(6));
t79 = t86 ^ 2;
t88 = cos(pkin(6));
t81 = 0.1e1 / t88 ^ 2;
t90 = cos(qJ(1));
t84 = t90 ^ 2;
t77 = t84 * t79 * t81 + 0.1e1;
t89 = sin(qJ(1));
t83 = t89 ^ 2;
t109 = 0.1e1 / t77 ^ 2 * t83;
t113 = t109 * t81;
t104 = t90 * t86;
t76 = atan2(t104, t88);
t72 = sin(t76);
t73 = cos(t76);
t58 = t72 * t104 + t73 * t88;
t55 = 0.1e1 / t58;
t87 = cos(pkin(12));
t103 = t90 * t87;
t85 = sin(pkin(12));
t107 = t89 * t85;
t99 = t88 * t107 - t103;
t65 = 0.1e1 / t99;
t80 = 0.1e1 / t88;
t56 = 0.1e1 / t58 ^ 2;
t66 = 0.1e1 / t99 ^ 2;
t112 = t56 * t89;
t105 = t90 * t85;
t106 = t89 * t87;
t70 = t88 * t106 + t105;
t111 = t66 * t70;
t69 = -t88 * t105 - t106;
t110 = t69 * t70;
t108 = t79 * t80;
t102 = qJD(1) * t90;
t74 = 0.1e1 / t77;
t101 = (t74 - 0.1e1) * t86;
t100 = -0.2e1 * t80 * t113;
t68 = t88 * t103 - t107;
t51 = (-t73 * t74 * t90 * t108 + t72 * t101) * t89;
t78 = t86 * t79;
t67 = t65 * t66;
t64 = t70 ^ 2;
t63 = t69 * qJD(1);
t62 = t68 * qJD(1);
t61 = t64 * t66 + 0.1e1;
t57 = t55 * t56;
t54 = t83 * t79 * t56 + 0.1e1;
t50 = qJD(1) * t51;
t1 = [(-t74 * t80 * t86 + t78 * t100) * t102, 0, 0, 0, 0, 0; (0.2e1 * (t51 * t112 - t55 * t90) / t54 ^ 2 * (-t50 * t57 * t83 + t102 * t112) * t79 + ((0.2e1 * t51 * t57 * t89 - t56 * t90) * t50 + (-t89 * t55 + ((-t51 + (-t78 * t113 - t101) * t89 * t72) * t90 - (t84 * t79 ^ 2 * t100 + (-t109 + (0.2e1 * t83 - t84) * t74) * t108) * t89 * t73) * t56) * qJD(1)) / t54) * t86, 0, 0, 0, 0, 0; 0.2e1 * (t66 * t110 + t65 * t68) / t61 ^ 2 * (t64 * t67 * t63 + t62 * t111) + (-t69 * t62 * t66 + (-0.2e1 * t67 * t110 - t68 * t66) * t63 + (-t99 * t111 + t70 * t65) * qJD(1)) / t61, 0, 0, 0, 0, 0;];
JaD_rot  = t1;
