% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR8_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobia_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:39
% EndTime: 2019-02-26 22:07:39
% DurationCPUTime: 0.15s
% Computational Cost: add. (359->30), mult. (1033->74), div. (77->9), fcn. (1491->11), ass. (0->49)
t67 = cos(pkin(6));
t69 = sin(qJ(2));
t73 = cos(qJ(1));
t76 = t73 * t69;
t70 = sin(qJ(1));
t72 = cos(qJ(2));
t77 = t70 * t72;
t60 = t67 * t76 + t77;
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t66 = sin(pkin(6));
t79 = t66 * t73;
t49 = t60 * t68 + t71 * t79;
t82 = t66 * t68;
t57 = -t67 * t71 + t69 * t82;
t46 = atan2(-t49, t57);
t42 = sin(t46);
t43 = cos(t46);
t41 = -t42 * t49 + t43 * t57;
t40 = 0.1e1 / t41 ^ 2;
t75 = t73 * t72;
t78 = t70 * t69;
t62 = -t67 * t78 + t75;
t81 = t66 * t71;
t52 = t62 * t68 - t70 * t81;
t88 = t40 * t52;
t87 = t43 * t49;
t53 = t62 * t71 + t70 * t82;
t48 = 0.1e1 / t53 ^ 2;
t61 = -t67 * t77 - t76;
t86 = t48 * t61;
t55 = 0.1e1 / t57 ^ 2;
t85 = t49 * t55;
t84 = t52 ^ 2 * t40;
t83 = t61 ^ 2 * t48;
t80 = t66 * t72;
t51 = t60 * t71 - t68 * t79;
t74 = -t42 * t57 - t87;
t59 = -t67 * t75 + t78;
t58 = t67 * t68 + t69 * t81;
t54 = 0.1e1 / t57;
t47 = 0.1e1 / t53;
t45 = 0.1e1 / (0.1e1 + t83);
t44 = 0.1e1 / (t49 ^ 2 * t55 + 0.1e1);
t39 = 0.1e1 / t41;
t38 = 0.1e1 / (0.1e1 + t84);
t37 = (t54 * t59 + t80 * t85) * t68 * t44;
t36 = (-t51 * t54 + t58 * t85) * t44;
t1 = [-t52 * t54 * t44, t37, t36, 0, 0, 0; (-t49 * t39 - (-t42 + (t54 * t87 + t42) * t44) * t84) * t38 (t61 * t68 * t39 - ((t42 * t59 + t43 * t80) * t68 + t74 * t37) * t88) * t38 (t53 * t39 - (t36 * t74 - t42 * t51 + t43 * t58) * t88) * t38, 0, 0, 0; (t59 * t47 + t51 * t86) * t45 (-t47 * t62 - t71 * t83) * t45, t52 * t45 * t86, 0, 0, 0;];
Ja_rot  = t1;
