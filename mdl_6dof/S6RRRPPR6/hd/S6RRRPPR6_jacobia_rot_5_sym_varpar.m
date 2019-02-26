% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR6_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:06:33
% EndTime: 2019-02-26 22:06:33
% DurationCPUTime: 0.13s
% Computational Cost: add. (727->32), mult. (1022->76), div. (77->9), fcn. (1479->11), ass. (0->49)
t74 = cos(pkin(6));
t75 = sin(qJ(2));
t78 = cos(qJ(1));
t81 = t78 * t75;
t76 = sin(qJ(1));
t77 = cos(qJ(2));
t82 = t76 * t77;
t65 = t74 * t81 + t82;
t72 = qJ(3) + pkin(11);
t70 = sin(t72);
t71 = cos(t72);
t73 = sin(pkin(6));
t84 = t73 * t78;
t53 = t65 * t70 + t71 * t84;
t87 = t73 * t75;
t60 = t70 * t87 - t74 * t71;
t51 = atan2(-t53, t60);
t48 = sin(t51);
t49 = cos(t51);
t47 = -t48 * t53 + t49 * t60;
t46 = 0.1e1 / t47 ^ 2;
t80 = t78 * t77;
t83 = t76 * t75;
t67 = -t74 * t83 + t80;
t86 = t73 * t76;
t56 = t67 * t70 - t71 * t86;
t92 = t46 * t56;
t91 = t49 * t53;
t59 = 0.1e1 / t60 ^ 2;
t90 = t53 * t59;
t89 = t56 ^ 2 * t46;
t57 = t67 * t71 + t70 * t86;
t66 = t74 * t82 + t81;
t63 = 0.1e1 / t66 ^ 2;
t88 = t57 * t63;
t85 = t73 * t77;
t55 = t65 * t71 - t70 * t84;
t79 = -t48 * t60 - t91;
t64 = t74 * t80 - t83;
t62 = 0.1e1 / t66;
t61 = t74 * t70 + t71 * t87;
t58 = 0.1e1 / t60;
t52 = 0.1e1 / (t57 ^ 2 * t63 + 0.1e1);
t50 = 0.1e1 / (t53 ^ 2 * t59 + 0.1e1);
t45 = 0.1e1 / t47;
t44 = 0.1e1 / (0.1e1 + t89);
t43 = (-t58 * t64 + t85 * t90) * t70 * t50;
t42 = (-t55 * t58 + t61 * t90) * t50;
t1 = [-t56 * t58 * t50, t43, t42, 0, 0, 0; (-t53 * t45 - (-t48 + (t58 * t91 + t48) * t50) * t89) * t44 (-t66 * t70 * t45 - ((-t48 * t64 + t49 * t85) * t70 + t79 * t43) * t92) * t44 (t57 * t45 - (t42 * t79 - t48 * t55 + t49 * t61) * t92) * t44, 0, 0, 0; (-t55 * t62 - t64 * t88) * t52 (-t62 * t66 * t71 - t67 * t88) * t52, -t56 * t62 * t52, 0, 0, 0;];
Ja_rot  = t1;
