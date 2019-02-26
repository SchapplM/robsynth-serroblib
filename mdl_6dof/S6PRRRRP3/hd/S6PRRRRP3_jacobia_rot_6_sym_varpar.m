% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRP3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:20
% EndTime: 2019-02-26 20:16:21
% DurationCPUTime: 0.14s
% Computational Cost: add. (427->31), mult. (1032->79), div. (70->9), fcn. (1461->13), ass. (0->51)
t88 = sin(pkin(6));
t93 = cos(qJ(3));
t100 = t88 * t93;
t87 = sin(pkin(11));
t89 = cos(pkin(11));
t94 = cos(qJ(2));
t90 = cos(pkin(6));
t92 = sin(qJ(2));
t98 = t90 * t92;
t78 = t87 * t94 + t89 * t98;
t91 = sin(qJ(3));
t70 = t100 * t89 + t78 * t91;
t101 = t88 * t91;
t81 = t101 * t92 - t90 * t93;
t69 = atan2(-t70, t81);
t66 = sin(t69);
t67 = cos(t69);
t60 = -t66 * t70 + t67 * t81;
t59 = 0.1e1 / t60 ^ 2;
t80 = -t87 * t98 + t89 * t94;
t73 = -t100 * t87 + t80 * t91;
t105 = t59 * t73;
t74 = t101 * t87 + t80 * t93;
t97 = t90 * t94;
t79 = t87 * t97 + t89 * t92;
t86 = qJ(4) + qJ(5);
t84 = sin(t86);
t85 = cos(t86);
t65 = t74 * t85 + t79 * t84;
t63 = 0.1e1 / t65 ^ 2;
t64 = t74 * t84 - t79 * t85;
t104 = t63 * t64;
t76 = 0.1e1 / t81 ^ 2;
t103 = t70 * t76;
t102 = t79 * t93;
t99 = t88 * t94;
t96 = t64 ^ 2 * t63 + 0.1e1;
t95 = -t66 * t81 - t67 * t70;
t82 = t100 * t92 + t90 * t91;
t77 = -t87 * t92 + t89 * t97;
t75 = 0.1e1 / t81;
t72 = -t101 * t89 + t78 * t93;
t68 = 0.1e1 / (t70 ^ 2 * t76 + 0.1e1);
t62 = 0.1e1 / t65;
t61 = 0.1e1 / t96;
t58 = 0.1e1 / t60;
t57 = 0.1e1 / (t73 ^ 2 * t59 + 0.1e1);
t56 = (t103 * t99 - t75 * t77) * t91 * t68;
t55 = (t103 * t82 - t72 * t75) * t68;
t54 = t96 * t61;
t1 = [0, t56, t55, 0, 0, 0; 0 (-t79 * t91 * t58 - ((-t66 * t77 + t67 * t99) * t91 + t95 * t56) * t105) * t57 (t74 * t58 - (t55 * t95 - t66 * t72 + t67 * t82) * t105) * t57, 0, 0, 0; 0 ((-t102 * t84 - t80 * t85) * t62 - (-t102 * t85 + t80 * t84) * t104) * t61 (t104 * t85 - t62 * t84) * t73 * t61, t54, t54, 0;];
Ja_rot  = t1;
