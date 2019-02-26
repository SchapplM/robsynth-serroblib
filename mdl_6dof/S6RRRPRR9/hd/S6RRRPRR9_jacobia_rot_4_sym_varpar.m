% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR9_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobia_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:27
% EndTime: 2019-02-26 22:20:27
% DurationCPUTime: 0.17s
% Computational Cost: add. (391->34), mult. (1119->83), div. (56->9), fcn. (1577->15), ass. (0->53)
t83 = cos(pkin(6));
t88 = cos(qJ(2));
t89 = cos(qJ(1));
t91 = t88 * t89;
t85 = sin(qJ(2));
t86 = sin(qJ(1));
t93 = t86 * t85;
t69 = -t83 * t91 + t93;
t79 = sin(pkin(7));
t82 = cos(pkin(7));
t80 = sin(pkin(6));
t95 = t80 * t89;
t59 = -t69 * t79 + t82 * t95;
t68 = -t79 * t80 * t88 + t82 * t83;
t58 = atan2(t59, t68);
t55 = sin(t58);
t56 = cos(t58);
t49 = t55 * t59 + t56 * t68;
t48 = 0.1e1 / t49 ^ 2;
t92 = t86 * t88;
t94 = t85 * t89;
t71 = -t83 * t92 - t94;
t96 = t80 * t86;
t60 = t71 * t79 - t82 * t96;
t101 = t48 * t60 ^ 2;
t78 = sin(pkin(13));
t81 = cos(pkin(13));
t84 = sin(qJ(3));
t87 = cos(qJ(3));
t90 = t78 * t87 + t81 * t84;
t63 = t90 * t79;
t65 = t90 * t82;
t72 = -t83 * t93 + t91;
t73 = t78 * t84 - t87 * t81;
t54 = t63 * t96 + t65 * t71 - t72 * t73;
t51 = 0.1e1 / t54 ^ 2;
t62 = t73 * t79;
t64 = t73 * t82;
t52 = -t62 * t96 - t64 * t71 - t72 * t90;
t100 = t51 * t52;
t99 = t52 ^ 2 * t51;
t98 = t56 * t59;
t97 = t80 * t85;
t70 = -t83 * t94 - t92;
t67 = 0.1e1 / t68 ^ 2;
t66 = 0.1e1 / t68;
t57 = 0.1e1 / (t59 ^ 2 * t67 + 0.1e1);
t50 = 0.1e1 / t54;
t47 = 0.1e1 / t49;
t46 = 0.1e1 / (0.1e1 + t101);
t45 = 0.1e1 / (0.1e1 + t99);
t44 = (-t59 * t67 * t97 + t66 * t70) * t79 * t57;
t1 = [t60 * t66 * t57, t44, 0, 0, 0, 0; (t59 * t47 + (t55 + (t66 * t98 - t55) * t57) * t101) * t46 (t72 * t79 * t47 + ((t55 * t70 + t56 * t97) * t79 + (-t55 * t68 + t98) * t44) * t60 * t48) * t46, 0, 0, 0, 0; ((t62 * t95 + t69 * t64 + t70 * t90) * t50 + (t63 * t95 + t69 * t65 - t70 * t73) * t100) * t45 ((-t64 * t72 + t71 * t90) * t50 + (-t65 * t72 - t71 * t73) * t100) * t45 (t54 * t50 + t99) * t45, 0, 0, 0;];
Ja_rot  = t1;
