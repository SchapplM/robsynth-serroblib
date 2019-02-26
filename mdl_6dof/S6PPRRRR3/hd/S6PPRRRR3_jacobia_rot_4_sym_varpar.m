% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPRRRR3
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRRR3_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobia_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:05
% EndTime: 2019-02-26 19:44:05
% DurationCPUTime: 0.16s
% Computational Cost: add. (558->35), mult. (1619->83), div. (35->9), fcn. (2209->17), ass. (0->50)
t70 = sin(pkin(14));
t75 = cos(pkin(14));
t76 = cos(pkin(13));
t71 = sin(pkin(13));
t79 = cos(pkin(6));
t93 = t71 * t79;
t69 = -t70 * t93 + t76 * t75;
t81 = sin(qJ(3));
t83 = cos(qJ(3));
t68 = -t76 * t70 - t75 * t93;
t78 = cos(pkin(7));
t73 = sin(pkin(7));
t74 = sin(pkin(6));
t92 = t73 * t74;
t84 = t68 * t78 + t71 * t92;
t61 = t69 * t83 + t84 * t81;
t80 = sin(qJ(4));
t95 = t61 * t80;
t82 = cos(qJ(4));
t94 = t61 * t82;
t91 = t73 * t79;
t90 = t74 * t78;
t89 = t75 * t78;
t88 = t76 * t79;
t60 = -t69 * t81 + t84 * t83;
t65 = -t68 * t73 + t71 * t90;
t72 = sin(pkin(8));
t77 = cos(pkin(8));
t86 = t60 * t77 + t65 * t72;
t51 = t86 * t80 + t94;
t49 = 0.1e1 / t51 ^ 2;
t50 = -t86 * t82 + t95;
t87 = t50 ^ 2 * t49 + 0.1e1;
t66 = -t71 * t70 + t75 * t88;
t85 = -t66 * t78 + t76 * t92;
t67 = t70 * t88 + t71 * t75;
t64 = -t81 * t91 + (-t70 * t83 - t81 * t89) * t74;
t62 = -(t83 * t91 + (-t70 * t81 + t83 * t89) * t74) * t72 + (-t75 * t92 + t79 * t78) * t77;
t59 = -t67 * t83 + t85 * t81;
t58 = 0.1e1 / t62 ^ 2;
t57 = -t60 * t72 + t65 * t77;
t56 = (-t67 * t81 - t85 * t83) * t72 - (-t66 * t73 - t76 * t90) * t77;
t55 = atan2(t56, t62);
t53 = cos(t55);
t52 = sin(t55);
t48 = 0.1e1 / t87;
t47 = t52 * t56 + t53 * t62;
t46 = 0.1e1 / t47 ^ 2;
t44 = (t59 / t62 + t64 * t56 * t58) * t72 / (t56 ^ 2 * t58 + 0.1e1);
t1 = [0, 0, t44, 0, 0, 0; 0, 0 (t61 * t72 / t47 - ((t52 * t59 - t53 * t64) * t72 + (-t52 * t62 + t53 * t56) * t44) * t57 * t46) / (t57 ^ 2 * t46 + 0.1e1) 0, 0, 0; 0, 0 ((t60 * t80 + t77 * t94) / t51 - (t60 * t82 - t77 * t95) * t50 * t49) * t48, t87 * t48, 0, 0;];
Ja_rot  = t1;
