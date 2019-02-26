% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPPRR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:35
% EndTime: 2019-02-26 19:44:35
% DurationCPUTime: 0.23s
% Computational Cost: add. (939->33), mult. (1727->83), div. (65->9), fcn. (2425->15), ass. (0->56)
t91 = sin(pkin(6));
t93 = cos(pkin(10));
t104 = t91 * t93;
t89 = sin(pkin(11));
t92 = cos(pkin(11));
t96 = sin(qJ(2));
t98 = cos(qJ(2));
t100 = t98 * t89 + t96 * t92;
t94 = cos(pkin(6));
t81 = t100 * t94;
t82 = t96 * t89 - t98 * t92;
t90 = sin(pkin(10));
t70 = t93 * t81 - t90 * t82;
t88 = pkin(12) + qJ(5);
t86 = sin(t88);
t87 = cos(t88);
t64 = t87 * t104 + t70 * t86;
t80 = t100 * t91;
t76 = t80 * t86 - t94 * t87;
t63 = atan2(-t64, t76);
t60 = sin(t63);
t61 = cos(t63);
t54 = -t60 * t64 + t61 * t76;
t53 = 0.1e1 / t54 ^ 2;
t101 = -t90 * t81 - t93 * t82;
t105 = t90 * t91;
t67 = t101 * t86 - t87 * t105;
t110 = t53 * t67;
t99 = t82 * t94;
t72 = -t100 * t93 + t90 * t99;
t95 = sin(qJ(6));
t107 = t72 * t95;
t68 = t101 * t87 + t86 * t105;
t97 = cos(qJ(6));
t59 = t68 * t97 - t107;
t57 = 0.1e1 / t59 ^ 2;
t106 = t72 * t97;
t58 = t68 * t95 + t106;
t109 = t57 * t58;
t75 = 0.1e1 / t76 ^ 2;
t108 = t64 * t75;
t103 = t58 ^ 2 * t57 + 0.1e1;
t102 = -t60 * t76 - t61 * t64;
t79 = t82 * t91;
t77 = t80 * t87 + t94 * t86;
t74 = 0.1e1 / t76;
t69 = -t100 * t90 - t93 * t99;
t66 = -t86 * t104 + t70 * t87;
t62 = 0.1e1 / (t64 ^ 2 * t75 + 0.1e1);
t56 = 0.1e1 / t59;
t55 = 0.1e1 / t103;
t52 = 0.1e1 / t54;
t51 = 0.1e1 / (t67 ^ 2 * t53 + 0.1e1);
t50 = (-t79 * t108 - t69 * t74) * t86 * t62;
t49 = (t77 * t108 - t66 * t74) * t62;
t1 = [0, t50, 0, 0, t49, 0; 0 (t72 * t86 * t52 - ((-t60 * t69 - t61 * t79) * t86 + t102 * t50) * t110) * t51, 0, 0 (t68 * t52 - (t102 * t49 - t60 * t66 + t61 * t77) * t110) * t51, 0; 0 ((-t101 * t97 + t87 * t107) * t56 - (t101 * t95 + t87 * t106) * t109) * t55, 0, 0 (t97 * t109 - t56 * t95) * t67 * t55, t103 * t55;];
Ja_rot  = t1;
