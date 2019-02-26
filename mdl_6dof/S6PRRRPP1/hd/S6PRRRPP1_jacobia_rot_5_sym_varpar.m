% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPP1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:44
% EndTime: 2019-02-26 20:08:44
% DurationCPUTime: 0.12s
% Computational Cost: add. (382->31), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->50)
t72 = sin(pkin(10));
t74 = cos(pkin(10));
t79 = cos(qJ(2));
t75 = cos(pkin(6));
t77 = sin(qJ(2));
t83 = t75 * t77;
t63 = t72 * t79 + t74 * t83;
t76 = sin(qJ(3));
t73 = sin(pkin(6));
t78 = cos(qJ(3));
t85 = t73 * t78;
t55 = t63 * t76 + t74 * t85;
t86 = t73 * t76;
t66 = -t75 * t78 + t77 * t86;
t54 = atan2(-t55, t66);
t51 = sin(t54);
t52 = cos(t54);
t45 = -t51 * t55 + t52 * t66;
t44 = 0.1e1 / t45 ^ 2;
t65 = -t72 * t83 + t74 * t79;
t58 = t65 * t76 - t72 * t85;
t90 = t44 * t58;
t59 = t65 * t78 + t72 * t86;
t82 = t75 * t79;
t64 = t72 * t82 + t74 * t77;
t71 = qJ(4) + pkin(11);
t69 = sin(t71);
t70 = cos(t71);
t50 = t59 * t70 + t64 * t69;
t48 = 0.1e1 / t50 ^ 2;
t49 = t59 * t69 - t64 * t70;
t89 = t48 * t49;
t61 = 0.1e1 / t66 ^ 2;
t88 = t55 * t61;
t87 = t64 * t78;
t84 = t73 * t79;
t81 = t49 ^ 2 * t48 + 0.1e1;
t80 = -t51 * t66 - t52 * t55;
t67 = t75 * t76 + t77 * t85;
t62 = -t72 * t77 + t74 * t82;
t60 = 0.1e1 / t66;
t57 = t63 * t78 - t74 * t86;
t53 = 0.1e1 / (t55 ^ 2 * t61 + 0.1e1);
t47 = 0.1e1 / t50;
t46 = 0.1e1 / t81;
t43 = 0.1e1 / t45;
t42 = 0.1e1 / (t58 ^ 2 * t44 + 0.1e1);
t41 = (-t60 * t62 + t84 * t88) * t76 * t53;
t40 = (-t57 * t60 + t67 * t88) * t53;
t1 = [0, t41, t40, 0, 0, 0; 0 (-t64 * t76 * t43 - ((-t51 * t62 + t52 * t84) * t76 + t80 * t41) * t90) * t42 (t59 * t43 - (t80 * t40 - t51 * t57 + t52 * t67) * t90) * t42, 0, 0, 0; 0 ((-t65 * t70 - t69 * t87) * t47 - (t65 * t69 - t70 * t87) * t89) * t46 (-t47 * t69 + t70 * t89) * t58 * t46, t81 * t46, 0, 0;];
Ja_rot  = t1;
