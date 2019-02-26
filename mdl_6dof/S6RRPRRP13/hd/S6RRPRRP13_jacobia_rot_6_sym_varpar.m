% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP13
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP13_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:52:49
% EndTime: 2019-02-26 21:52:49
% DurationCPUTime: 0.15s
% Computational Cost: add. (459->36), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->55)
t78 = cos(pkin(6));
t85 = cos(qJ(2));
t86 = cos(qJ(1));
t91 = t86 * t85;
t81 = sin(qJ(2));
t82 = sin(qJ(1));
t94 = t82 * t81;
t73 = -t78 * t91 + t94;
t80 = sin(qJ(4));
t84 = cos(qJ(4));
t77 = sin(pkin(6));
t95 = t77 * t86;
t65 = t73 * t84 + t80 * t95;
t96 = t77 * t85;
t71 = t78 * t80 + t84 * t96;
t62 = atan2(t65, t71);
t59 = sin(t62);
t60 = cos(t62);
t53 = t59 * t65 + t60 * t71;
t52 = 0.1e1 / t53 ^ 2;
t92 = t86 * t81;
t93 = t82 * t85;
t87 = t78 * t93 + t92;
t97 = t77 * t82;
t63 = t80 * t97 - t87 * t84;
t105 = t52 * t63;
t75 = -t78 * t94 + t91;
t79 = sin(qJ(5));
t100 = t75 * t79;
t64 = t87 * t80 + t84 * t97;
t83 = cos(qJ(5));
t58 = t64 * t83 + t100;
t56 = 0.1e1 / t58 ^ 2;
t99 = t75 * t83;
t57 = t64 * t79 - t99;
t104 = t56 * t57;
t103 = t60 * t65;
t102 = t63 ^ 2 * t52;
t70 = 0.1e1 / t71 ^ 2;
t101 = t65 * t70;
t98 = t77 * t81;
t90 = t57 ^ 2 * t56 + 0.1e1;
t89 = -t59 * t71 + t103;
t88 = -t73 * t80 + t84 * t95;
t74 = t78 * t92 + t93;
t72 = t78 * t84 - t80 * t96;
t69 = 0.1e1 / t71;
t61 = 0.1e1 / (t65 ^ 2 * t70 + 0.1e1);
t55 = 0.1e1 / t58;
t54 = 0.1e1 / t90;
t51 = 0.1e1 / t53;
t50 = 0.1e1 / (0.1e1 + t102);
t49 = (t98 * t101 + t69 * t74) * t84 * t61;
t48 = (-t72 * t101 + t69 * t88) * t61;
t1 = [-t63 * t69 * t61, t49, 0, t48, 0, 0; (t65 * t51 - (-t59 + (-t69 * t103 + t59) * t61) * t102) * t50 (-t75 * t84 * t51 - ((t59 * t74 - t60 * t98) * t84 + t89 * t49) * t105) * t50, 0 (t64 * t51 - (t89 * t48 + t59 * t88 + t60 * t72) * t105) * t50, 0, 0; ((t74 * t83 + t79 * t88) * t55 - (-t74 * t79 + t83 * t88) * t104) * t54 ((t80 * t100 + t87 * t83) * t55 - (-t87 * t79 + t80 * t99) * t104) * t54, 0 (t83 * t104 - t79 * t55) * t63 * t54, t90 * t54, 0;];
Ja_rot  = t1;
