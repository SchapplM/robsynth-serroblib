% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:13:09
% EndTime: 2019-02-26 22:13:09
% DurationCPUTime: 0.10s
% Computational Cost: add. (141->25), mult. (425->68), div. (58->9), fcn. (611->11), ass. (0->40)
t70 = sin(qJ(1));
t72 = cos(qJ(3));
t73 = cos(qJ(2));
t68 = sin(qJ(3));
t74 = cos(qJ(1));
t77 = t74 * t68;
t57 = -t70 * t72 + t73 * t77;
t76 = t74 * t72;
t58 = t70 * t68 + t73 * t76;
t67 = sin(qJ(5));
t71 = cos(qJ(5));
t50 = t57 * t67 + t58 * t71;
t48 = 0.1e1 / t50 ^ 2;
t49 = -t57 * t71 + t58 * t67;
t85 = t48 * t49;
t84 = t49 ^ 2 * t48;
t69 = sin(qJ(2));
t79 = t70 * t69;
t62 = atan2(t79, t73);
t59 = sin(t62);
t60 = cos(t62);
t53 = t59 * t79 + t60 * t73;
t52 = 0.1e1 / t53 ^ 2;
t83 = t52 * t74 ^ 2;
t64 = t69 ^ 2;
t82 = t64 / t73 ^ 2;
t81 = t69 * t74;
t61 = 0.1e1 / (t70 ^ 2 * t82 + 0.1e1);
t80 = t70 * t61;
t78 = t70 * t73;
t75 = 0.1e1 + t84;
t65 = 0.1e1 / t73;
t56 = -t72 * t78 + t77;
t55 = -t68 * t78 - t76;
t54 = (0.1e1 + t82) * t80;
t51 = 0.1e1 / t53;
t47 = 0.1e1 / t50;
t46 = 0.1e1 / (t64 * t83 + 0.1e1);
t45 = 0.1e1 / t75;
t1 = [t65 * t61 * t81, t54, 0, 0, 0, 0; (t51 * t79 + (t60 * t64 * t65 * t80 + (-t61 + 0.1e1) * t69 * t59) * t69 * t83) * t46 (-t73 * t51 + (t59 * t78 - t60 * t69 + (-t59 * t73 + t60 * t79) * t54) * t69 * t52) * t74 * t46, 0, 0, 0, 0; ((-t55 * t71 + t56 * t67) * t47 - (t55 * t67 + t56 * t71) * t85) * t45 ((-t67 * t72 + t68 * t71) * t47 - (-t67 * t68 - t71 * t72) * t85) * t45 * t81 (-t50 * t47 - t84) * t45, 0, t75 * t45, 0;];
Ja_rot  = t1;
