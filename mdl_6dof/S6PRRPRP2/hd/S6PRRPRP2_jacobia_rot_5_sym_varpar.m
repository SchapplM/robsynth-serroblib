% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRP2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:47
% EndTime: 2019-02-26 20:01:47
% DurationCPUTime: 0.18s
% Computational Cost: add. (641->31), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->52)
t71 = sin(pkin(10));
t73 = cos(pkin(10));
t78 = cos(qJ(2));
t74 = cos(pkin(6));
t76 = sin(qJ(2));
t82 = t74 * t76;
t64 = t71 * t78 + t73 * t82;
t70 = qJ(3) + pkin(11);
t68 = sin(t70);
t69 = cos(t70);
t72 = sin(pkin(6));
t85 = t72 * t73;
t54 = t64 * t68 + t69 * t85;
t84 = t72 * t76;
t61 = t68 * t84 - t69 * t74;
t53 = atan2(-t54, t61);
t48 = sin(t53);
t49 = cos(t53);
t44 = -t48 * t54 + t49 * t61;
t43 = 0.1e1 / t44 ^ 2;
t66 = -t71 * t82 + t73 * t78;
t86 = t71 * t72;
t57 = t66 * t68 - t69 * t86;
t91 = t43 * t57;
t58 = t66 * t69 + t68 * t86;
t77 = cos(qJ(5));
t81 = t74 * t78;
t65 = t71 * t81 + t73 * t76;
t75 = sin(qJ(5));
t88 = t65 * t75;
t51 = t58 * t77 + t88;
t47 = 0.1e1 / t51 ^ 2;
t87 = t65 * t77;
t50 = t58 * t75 - t87;
t90 = t47 * t50;
t60 = 0.1e1 / t61 ^ 2;
t89 = t54 * t60;
t83 = t72 * t78;
t80 = t47 * t50 ^ 2 + 0.1e1;
t79 = -t48 * t61 - t49 * t54;
t63 = -t71 * t76 + t73 * t81;
t62 = t68 * t74 + t69 * t84;
t59 = 0.1e1 / t61;
t56 = t64 * t69 - t68 * t85;
t52 = 0.1e1 / (t54 ^ 2 * t60 + 0.1e1);
t46 = 0.1e1 / t51;
t45 = 0.1e1 / t80;
t42 = 0.1e1 / t44;
t41 = 0.1e1 / (t43 * t57 ^ 2 + 0.1e1);
t40 = (-t59 * t63 + t83 * t89) * t68 * t52;
t39 = (-t56 * t59 + t62 * t89) * t52;
t1 = [0, t40, t39, 0, 0, 0; 0 (-t65 * t68 * t42 - ((-t48 * t63 + t49 * t83) * t68 + t79 * t40) * t91) * t41 (t58 * t42 - (t39 * t79 - t48 * t56 + t49 * t62) * t91) * t41, 0, 0, 0; 0 ((-t66 * t77 - t69 * t88) * t46 - (t66 * t75 - t69 * t87) * t90) * t45 (-t46 * t75 + t77 * t90) * t57 * t45, 0, t80 * t45, 0;];
Ja_rot  = t1;
