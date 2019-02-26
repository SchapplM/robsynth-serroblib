% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP14
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

function Ja_rot = S6RRPRRP14_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:53:34
% EndTime: 2019-02-26 21:53:34
% DurationCPUTime: 0.16s
% Computational Cost: add. (459->36), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->55)
t72 = cos(pkin(6));
t79 = cos(qJ(2));
t80 = cos(qJ(1));
t85 = t80 * t79;
t75 = sin(qJ(2));
t76 = sin(qJ(1));
t88 = t76 * t75;
t67 = -t72 * t85 + t88;
t74 = sin(qJ(4));
t78 = cos(qJ(4));
t71 = sin(pkin(6));
t89 = t71 * t80;
t59 = t67 * t78 + t74 * t89;
t90 = t71 * t79;
t65 = t72 * t74 + t78 * t90;
t56 = atan2(t59, t65);
t53 = sin(t56);
t54 = cos(t56);
t47 = t53 * t59 + t54 * t65;
t46 = 0.1e1 / t47 ^ 2;
t86 = t80 * t75;
t87 = t76 * t79;
t81 = t72 * t87 + t86;
t91 = t71 * t76;
t57 = t74 * t91 - t81 * t78;
t99 = t46 * t57;
t58 = t81 * t74 + t78 * t91;
t77 = cos(qJ(5));
t69 = -t72 * t88 + t85;
t73 = sin(qJ(5));
t94 = t69 * t73;
t52 = t58 * t77 + t94;
t50 = 0.1e1 / t52 ^ 2;
t93 = t69 * t77;
t51 = t58 * t73 - t93;
t98 = t50 * t51;
t97 = t54 * t59;
t96 = t57 ^ 2 * t46;
t64 = 0.1e1 / t65 ^ 2;
t95 = t59 * t64;
t92 = t71 * t75;
t84 = t51 ^ 2 * t50 + 0.1e1;
t83 = -t53 * t65 + t97;
t82 = -t67 * t74 + t78 * t89;
t68 = t72 * t86 + t87;
t66 = t72 * t78 - t74 * t90;
t63 = 0.1e1 / t65;
t55 = 0.1e1 / (t59 ^ 2 * t64 + 0.1e1);
t49 = 0.1e1 / t52;
t48 = 0.1e1 / t84;
t45 = 0.1e1 / t47;
t44 = 0.1e1 / (0.1e1 + t96);
t43 = (t63 * t68 + t92 * t95) * t78 * t55;
t42 = (t63 * t82 - t66 * t95) * t55;
t1 = [-t57 * t63 * t55, t43, 0, t42, 0, 0; (t59 * t45 - (-t53 + (-t63 * t97 + t53) * t55) * t96) * t44 (-t69 * t78 * t45 - ((t53 * t68 - t54 * t92) * t78 + t83 * t43) * t99) * t44, 0 (t58 * t45 - (t83 * t42 + t53 * t82 + t54 * t66) * t99) * t44, 0, 0; ((t68 * t77 + t73 * t82) * t49 - (-t68 * t73 + t77 * t82) * t98) * t48 ((t74 * t94 + t81 * t77) * t49 - (-t81 * t73 + t74 * t93) * t98) * t48, 0 (-t73 * t49 + t77 * t98) * t57 * t48, t84 * t48, 0;];
Ja_rot  = t1;
