% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR11_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:49
% EndTime: 2019-02-26 22:21:50
% DurationCPUTime: 0.17s
% Computational Cost: add. (270->29), mult. (779->75), div. (78->11), fcn. (1122->13), ass. (0->49)
t70 = cos(pkin(6));
t77 = cos(qJ(2));
t78 = cos(qJ(1));
t80 = t78 * t77;
t73 = sin(qJ(2));
t74 = sin(qJ(1));
t83 = t74 * t73;
t63 = -t70 * t83 + t80;
t72 = sin(qJ(3));
t76 = cos(qJ(3));
t69 = sin(pkin(6));
t86 = t69 * t74;
t53 = t63 * t72 - t76 * t86;
t54 = t63 * t76 + t72 * t86;
t71 = sin(qJ(5));
t75 = cos(qJ(5));
t50 = t53 * t71 + t54 * t75;
t48 = 0.1e1 / t50 ^ 2;
t49 = -t53 * t75 + t54 * t71;
t91 = t48 * t49;
t90 = t49 ^ 2 * t48;
t59 = -t70 * t80 + t83;
t85 = t69 * t77;
t58 = atan2(t59, t85);
t56 = cos(t58);
t89 = t56 * t59;
t55 = sin(t58);
t46 = t55 * t59 + t56 * t85;
t45 = 0.1e1 / t46 ^ 2;
t81 = t78 * t73;
t82 = t74 * t77;
t61 = t70 * t82 + t81;
t88 = t61 ^ 2 * t45;
t66 = 0.1e1 / t69;
t67 = 0.1e1 / t77;
t87 = t66 * t67;
t84 = t69 * t78;
t79 = 0.1e1 + t90;
t68 = 0.1e1 / t77 ^ 2;
t60 = t70 * t81 + t82;
t57 = 0.1e1 / (0.1e1 + t59 ^ 2 / t69 ^ 2 * t68);
t52 = -t60 * t76 + t72 * t84;
t51 = -t60 * t72 - t76 * t84;
t47 = 0.1e1 / t50;
t44 = 0.1e1 / t46;
t43 = 0.1e1 / (0.1e1 + t88);
t42 = (t59 * t68 * t73 + t60 * t67) * t66 * t57;
t41 = 0.1e1 / t79;
t1 = [t61 * t57 * t87, t42, 0, 0, 0, 0; (t59 * t44 + (t55 + (t87 * t89 - t55) * t57) * t88) * t43 (-t63 * t44 + (-t56 * t69 * t73 + t55 * t60 + (-t55 * t85 + t89) * t42) * t61 * t45) * t43, 0, 0, 0, 0; ((-t51 * t75 + t52 * t71) * t47 - (t51 * t71 + t52 * t75) * t91) * t41 ((-t71 * t76 + t72 * t75) * t47 - (-t71 * t72 - t75 * t76) * t91) * t41 * t61 (-t50 * t47 - t90) * t41, 0, t79 * t41, 0;];
Ja_rot  = t1;
