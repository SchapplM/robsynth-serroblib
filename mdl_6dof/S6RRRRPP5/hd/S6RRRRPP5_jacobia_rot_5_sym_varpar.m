% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:27:43
% EndTime: 2019-02-26 22:27:43
% DurationCPUTime: 0.14s
% Computational Cost: add. (610->28), mult. (689->68), div. (139->11), fcn. (1046->9), ass. (0->42)
t75 = sin(qJ(2));
t90 = t75 ^ 2;
t74 = qJ(3) + qJ(4);
t68 = sin(t74);
t69 = cos(t74);
t78 = cos(qJ(1));
t80 = t78 * t69;
t76 = sin(qJ(1));
t77 = cos(qJ(2));
t82 = t76 * t77;
t59 = t68 * t82 + t80;
t84 = t75 * t68;
t55 = atan2(-t59, t84);
t52 = sin(t55);
t53 = cos(t55);
t50 = -t52 * t59 + t53 * t84;
t49 = 0.1e1 / t50 ^ 2;
t81 = t78 * t68;
t62 = -t76 * t69 + t77 * t81;
t89 = t49 * t62;
t87 = t53 * t59;
t86 = t62 ^ 2 * t49;
t66 = 0.1e1 / t68;
t71 = 0.1e1 / t75;
t85 = t66 * t71;
t83 = t75 * t78;
t63 = t76 * t68 + t77 * t80;
t58 = 0.1e1 / t63 ^ 2;
t79 = t78 ^ 2 * t90 * t58;
t72 = 0.1e1 / t90;
t67 = 0.1e1 / t68 ^ 2;
t61 = t69 * t82 - t81;
t57 = 0.1e1 / t63;
t56 = 0.1e1 / (0.1e1 + t79);
t54 = 0.1e1 / (t59 ^ 2 * t72 * t67 + 0.1e1);
t51 = t62 * t58 * t56 * t83;
t48 = 0.1e1 / t50;
t47 = (t59 * t66 * t72 * t77 + t76) * t54;
t46 = 0.1e1 / (0.1e1 + t86);
t45 = (t59 * t67 * t69 - t61 * t66) * t71 * t54;
t44 = (t63 * t48 - (t53 * t75 * t69 - t52 * t61 + (-t52 * t84 - t87) * t45) * t89) * t46;
t1 = [-t62 * t54 * t85, t47, t45, t45, 0, 0; (-t59 * t48 - (-t52 + (t85 * t87 + t52) * t54) * t86) * t46 (t47 * t87 * t89 + (-t48 * t83 - (t53 * t77 + (-t47 + t76) * t52 * t75) * t89) * t68) * t46, t44, t44, 0, 0; (-t58 * t61 * t78 + t57 * t76) * t75 * t56 (-t57 * t77 * t78 - t69 * t79) * t56, -t51, -t51, 0, 0;];
Ja_rot  = t1;
