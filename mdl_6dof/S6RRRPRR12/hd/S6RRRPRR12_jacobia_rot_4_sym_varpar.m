% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR12
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR12_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_jacobia_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:22:25
% EndTime: 2019-02-26 22:22:26
% DurationCPUTime: 0.18s
% Computational Cost: add. (428->37), mult. (1217->90), div. (80->9), fcn. (1746->13), ass. (0->53)
t72 = cos(pkin(6));
t74 = sin(qJ(2));
t78 = cos(qJ(1));
t81 = t78 * t74;
t75 = sin(qJ(1));
t77 = cos(qJ(2));
t82 = t75 * t77;
t64 = t72 * t81 + t82;
t73 = sin(qJ(3));
t76 = cos(qJ(3));
t70 = sin(pkin(6));
t84 = t70 * t78;
t53 = t64 * t73 + t76 * t84;
t87 = t70 * t73;
t61 = -t72 * t76 + t74 * t87;
t52 = atan2(-t53, t61);
t49 = sin(t52);
t50 = cos(t52);
t43 = -t49 * t53 + t50 * t61;
t42 = 0.1e1 / t43 ^ 2;
t80 = t78 * t77;
t83 = t75 * t74;
t66 = -t72 * t83 + t80;
t86 = t70 * t76;
t57 = t66 * t73 - t75 * t86;
t93 = t42 * t57;
t58 = t66 * t76 + t75 * t87;
t65 = t72 * t82 + t81;
t69 = sin(pkin(12));
t71 = cos(pkin(12));
t48 = t58 * t71 + t65 * t69;
t46 = 0.1e1 / t48 ^ 2;
t47 = t58 * t69 - t65 * t71;
t92 = t46 * t47;
t91 = t50 * t53;
t60 = 0.1e1 / t61 ^ 2;
t90 = t53 * t60;
t89 = t57 ^ 2 * t42;
t88 = t65 * t76;
t85 = t70 * t77;
t55 = t64 * t76 - t73 * t84;
t79 = -t49 * t61 - t91;
t63 = t72 * t80 - t83;
t62 = t72 * t73 + t74 * t86;
t59 = 0.1e1 / t61;
t51 = 0.1e1 / (t53 ^ 2 * t60 + 0.1e1);
t45 = 0.1e1 / t48;
t44 = 0.1e1 / (t47 ^ 2 * t46 + 0.1e1);
t41 = 0.1e1 / t43;
t40 = 0.1e1 / (0.1e1 + t89);
t39 = (-t59 * t63 + t85 * t90) * t73 * t51;
t38 = (-t55 * t59 + t62 * t90) * t51;
t1 = [-t57 * t59 * t51, t39, t38, 0, 0, 0; (-t53 * t41 - (-t49 + (t59 * t91 + t49) * t51) * t89) * t40 (-t65 * t73 * t41 - ((-t49 * t63 + t50 * t85) * t73 + t79 * t39) * t93) * t40 (t58 * t41 - (t79 * t38 - t49 * t55 + t50 * t62) * t93) * t40, 0, 0, 0; ((-t55 * t69 - t63 * t71) * t45 - (-t55 * t71 + t63 * t69) * t92) * t44 ((-t66 * t71 - t69 * t88) * t45 - (t66 * t69 - t71 * t88) * t92) * t44 (-t69 * t45 + t71 * t92) * t57 * t44, 0, 0, 0;];
Ja_rot  = t1;
