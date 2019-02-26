% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR10_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:47
% EndTime: 2019-02-26 22:08:47
% DurationCPUTime: 0.17s
% Computational Cost: add. (428->37), mult. (1217->90), div. (80->9), fcn. (1746->13), ass. (0->53)
t69 = cos(pkin(6));
t71 = sin(qJ(2));
t75 = cos(qJ(1));
t78 = t75 * t71;
t72 = sin(qJ(1));
t74 = cos(qJ(2));
t79 = t72 * t74;
t62 = t69 * t78 + t79;
t70 = sin(qJ(3));
t73 = cos(qJ(3));
t67 = sin(pkin(6));
t81 = t67 * t75;
t52 = t62 * t73 - t70 * t81;
t83 = t67 * t73;
t60 = t69 * t70 + t71 * t83;
t50 = atan2(-t52, t60);
t47 = sin(t50);
t48 = cos(t50);
t41 = -t47 * t52 + t48 * t60;
t40 = 0.1e1 / t41 ^ 2;
t77 = t75 * t74;
t80 = t72 * t71;
t64 = -t69 * t80 + t77;
t84 = t67 * t70;
t56 = t64 * t73 + t72 * t84;
t90 = t40 * t56;
t55 = t64 * t70 - t72 * t83;
t63 = t69 * t79 + t78;
t66 = sin(pkin(11));
t68 = cos(pkin(11));
t46 = t55 * t66 + t63 * t68;
t44 = 0.1e1 / t46 ^ 2;
t45 = -t55 * t68 + t63 * t66;
t89 = t44 * t45;
t88 = t48 * t52;
t58 = 0.1e1 / t60 ^ 2;
t87 = t52 * t58;
t86 = t56 ^ 2 * t40;
t85 = t63 * t70;
t82 = t67 * t74;
t76 = -t47 * t60 - t88;
t51 = t62 * t70 + t73 * t81;
t61 = t69 * t77 - t80;
t59 = t69 * t73 - t71 * t84;
t57 = 0.1e1 / t60;
t49 = 0.1e1 / (t52 ^ 2 * t58 + 0.1e1);
t43 = 0.1e1 / t46;
t42 = 0.1e1 / (t45 ^ 2 * t44 + 0.1e1);
t39 = 0.1e1 / t41;
t38 = 0.1e1 / (0.1e1 + t86);
t37 = (-t57 * t61 + t82 * t87) * t73 * t49;
t36 = (t51 * t57 + t59 * t87) * t49;
t1 = [-t56 * t57 * t49, t37, t36, 0, 0, 0; (-t52 * t39 - (-t47 + (t57 * t88 + t47) * t49) * t86) * t38 (-t63 * t73 * t39 - ((-t47 * t61 + t48 * t82) * t73 + t76 * t37) * t90) * t38 (-t55 * t39 - (t76 * t36 + t47 * t51 + t48 * t59) * t90) * t38, 0, 0, 0; ((t51 * t68 + t61 * t66) * t43 - (-t51 * t66 + t61 * t68) * t89) * t42 ((t64 * t66 + t68 * t85) * t43 - (t64 * t68 - t66 * t85) * t89) * t42 (-t68 * t43 - t66 * t89) * t56 * t42, 0, 0, 0;];
Ja_rot  = t1;
