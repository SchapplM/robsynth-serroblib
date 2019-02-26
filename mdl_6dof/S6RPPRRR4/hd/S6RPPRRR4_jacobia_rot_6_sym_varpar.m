% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRR4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:36:28
% EndTime: 2019-02-26 20:36:28
% DurationCPUTime: 0.12s
% Computational Cost: add. (282->21), mult. (497->61), div. (57->9), fcn. (736->11), ass. (0->40)
t76 = sin(pkin(10));
t77 = cos(pkin(10));
t87 = sin(qJ(1));
t88 = cos(qJ(1));
t60 = -t87 * t76 - t77 * t88;
t89 = t60 ^ 2;
t74 = cos(qJ(4));
t61 = t76 * t88 - t87 * t77;
t73 = sin(qJ(4));
t81 = t61 * t73;
t58 = atan2(t81, t74);
t56 = sin(t58);
t57 = cos(t58);
t51 = t56 * t81 + t57 * t74;
t50 = 0.1e1 / t51 ^ 2;
t86 = t50 * t73;
t72 = qJ(5) + qJ(6);
t67 = sin(t72);
t68 = cos(t72);
t79 = t68 * t74;
t55 = -t60 * t79 + t61 * t67;
t53 = 0.1e1 / t55 ^ 2;
t80 = t67 * t74;
t54 = -t60 * t80 - t61 * t68;
t85 = t53 * t54;
t84 = t56 * t74;
t83 = t60 * t73;
t69 = t73 ^ 2;
t78 = t69 / t74 ^ 2;
t59 = 0.1e1 / (t61 ^ 2 * t78 + 0.1e1);
t82 = t61 * t59;
t75 = t54 ^ 2 * t53 + 0.1e1;
t70 = 0.1e1 / t74;
t52 = 0.1e1 / t55;
t49 = 0.1e1 / t51;
t48 = (0.1e1 + t78) * t82;
t47 = 0.1e1 / t75;
t46 = 0.1e1 / (t89 * t69 * t50 + 0.1e1);
t45 = t75 * t47;
t1 = [t70 * t59 * t83, 0, 0, t48, 0, 0; (t49 * t81 + (t57 * t69 * t70 * t82 + (-t59 + 0.1e1) * t73 * t56) * t89 * t86) * t46, 0, 0 (-t74 * t49 + (t61 * t84 - t57 * t73 + (t57 * t81 - t84) * t48) * t86) * t60 * t46, 0, 0; ((-t60 * t68 + t61 * t80) * t52 - (t60 * t67 + t61 * t79) * t85) * t47, 0, 0 (t52 * t67 - t68 * t85) * t47 * t83, t45, t45;];
Ja_rot  = t1;
