% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRPRRR5
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
% Datum: 2019-02-26 19:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR5_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobiaD_rot_3_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:11
% EndTime: 2019-02-26 19:56:11
% DurationCPUTime: 0.27s
% Computational Cost: add. (536->40), mult. (1574->99), div. (384->14), fcn. (2127->9), ass. (0->50)
t90 = sin(pkin(6));
t85 = 0.1e1 / t90 ^ 2;
t89 = sin(pkin(11));
t118 = 0.1e1 / t89 ^ 2 * t85;
t92 = cos(qJ(2));
t110 = t90 * t92;
t108 = cos(pkin(11));
t109 = cos(pkin(6));
t104 = t109 * t108;
t91 = sin(qJ(2));
t78 = -t104 * t92 + t89 * t91;
t67 = atan2(-t78, -t110);
t65 = sin(t67);
t66 = cos(t67);
t63 = -t110 * t66 - t65 * t78;
t60 = 0.1e1 / t63;
t84 = 0.1e1 / t90;
t86 = 0.1e1 / t92;
t61 = 0.1e1 / t63 ^ 2;
t87 = 0.1e1 / t92 ^ 2;
t105 = t89 * t109;
t101 = -t105 * t92 - t108 * t91;
t115 = -0.2e1 * t101;
t103 = t110 * t65 - t66 * t78;
t107 = t66 * t90 * t91;
t111 = t87 * t91;
t106 = t78 * t111;
t76 = t78 ^ 2;
t71 = t76 * t85 * t87 + 0.1e1;
t68 = 0.1e1 / t71;
t112 = t68 * t84;
t80 = t104 * t91 + t89 * t92;
t73 = t80 * qJD(2);
t55 = (qJD(2) * t106 + t73 * t86) * t112;
t53 = qJD(2) * t107 + t103 * t55 - t65 * t73;
t114 = t53 * t60 * t61;
t113 = t61 * t101;
t102 = t80 * t86 + t106;
t82 = -t105 * t91 + t108 * t92;
t88 = t86 * t87;
t77 = t101 ^ 2;
t75 = t82 * qJD(2);
t74 = t101 * qJD(2);
t72 = qJD(2) * t78;
t70 = t118 * t82 ^ 2 + 0.1e1;
t59 = t61 * t77 + 0.1e1;
t56 = t102 * t112;
t54 = t103 * t56 - t65 * t80 + t107;
t52 = (-0.2e1 * t102 / t71 ^ 2 * (qJD(2) * t76 * t88 * t91 + t73 * t78 * t87) * t85 + (t73 * t111 - t72 * t86 + (t80 * t111 + (0.2e1 * t88 * t91 ^ 2 + t86) * t78) * qJD(2)) * t68) * t84;
t1 = [0, t52, 0, 0, 0, 0; 0, 0.2e1 * (-t113 * t54 - t60 * t82) / t59 ^ 2 * (-t113 * t75 - t114 * t77) + (t54 * t114 * t115 + t74 * t60 + (-t82 * t53 - t54 * t75 - (-(qJD(2) * t110 - t52 * t78 - t56 * t73 + (t110 * t56 - t80) * t55) * t66 - (t55 * t56 * t78 + t72 + (t52 * t92 + (-qJD(2) * t56 - t55) * t91) * t90) * t65) * t101) * t61) / t59, 0, 0, 0, 0; 0 (-t75 / t70 + 0.1e1 / t70 ^ 2 * t82 * t74 * t115 * t118) * t84 / t89, 0, 0, 0, 0;];
JaD_rot  = t1;
