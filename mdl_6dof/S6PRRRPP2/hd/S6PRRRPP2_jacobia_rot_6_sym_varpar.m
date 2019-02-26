% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPP2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:09:33
% EndTime: 2019-02-26 20:09:33
% DurationCPUTime: 0.12s
% Computational Cost: add. (338->30), mult. (960->79), div. (66->9), fcn. (1365->13), ass. (0->49)
t81 = sin(pkin(10));
t83 = cos(pkin(10));
t90 = cos(qJ(2));
t84 = cos(pkin(6));
t87 = sin(qJ(2));
t93 = t84 * t87;
t76 = t81 * t90 + t83 * t93;
t86 = sin(qJ(3));
t82 = sin(pkin(6));
t89 = cos(qJ(3));
t95 = t82 * t89;
t69 = t76 * t86 + t83 * t95;
t96 = t82 * t86;
t79 = t84 * t89 - t87 * t96;
t68 = atan2(t69, t79);
t65 = sin(t68);
t66 = cos(t68);
t58 = t65 * t69 + t66 * t79;
t57 = 0.1e1 / t58 ^ 2;
t78 = -t81 * t93 + t83 * t90;
t71 = -t78 * t86 + t81 * t95;
t101 = t57 * t71;
t72 = t78 * t89 + t81 * t96;
t92 = t84 * t90;
t77 = t81 * t92 + t83 * t87;
t85 = sin(qJ(4));
t88 = cos(qJ(4));
t64 = t72 * t88 + t77 * t85;
t62 = 0.1e1 / t64 ^ 2;
t63 = -t72 * t85 + t77 * t88;
t100 = t63 ^ 2 * t62;
t99 = t62 * t63;
t74 = 0.1e1 / t79 ^ 2;
t98 = t69 * t74;
t97 = t77 * t89;
t94 = t82 * t90;
t91 = -t65 * t79 + t66 * t69;
t80 = -t84 * t86 - t87 * t95;
t75 = -t81 * t87 + t83 * t92;
t73 = 0.1e1 / t79;
t70 = t76 * t89 - t83 * t96;
t67 = 0.1e1 / (t69 ^ 2 * t74 + 0.1e1);
t61 = 0.1e1 / t64;
t59 = 0.1e1 / (0.1e1 + t100);
t56 = 0.1e1 / t58;
t55 = 0.1e1 / (t57 * t71 ^ 2 + 0.1e1);
t54 = (t73 * t75 + t94 * t98) * t86 * t67;
t53 = (t70 * t73 - t80 * t98) * t67;
t1 = [0, t54, t53, 0, 0, 0; 0 (t77 * t86 * t56 - ((t65 * t75 - t66 * t94) * t86 + t91 * t54) * t101) * t55 (-t72 * t56 - (t91 * t53 + t65 * t70 + t66 * t80) * t101) * t55, 0, 0, 0; 0 ((t78 * t88 + t85 * t97) * t61 - (t78 * t85 - t88 * t97) * t99) * t59 (-t61 * t85 - t88 * t99) * t71 * t59 (-t61 * t64 - t100) * t59, 0, 0;];
Ja_rot  = t1;
