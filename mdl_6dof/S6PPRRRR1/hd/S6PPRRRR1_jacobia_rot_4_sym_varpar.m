% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPRRRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRRR1_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobia_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:42
% EndTime: 2019-02-26 19:42:43
% DurationCPUTime: 0.14s
% Computational Cost: add. (333->29), mult. (995->67), div. (35->9), fcn. (1360->15), ass. (0->42)
t64 = sin(pkin(13));
t65 = sin(pkin(12));
t68 = cos(pkin(13));
t69 = cos(pkin(12));
t70 = cos(pkin(7));
t71 = cos(pkin(6));
t81 = t69 * t71;
t66 = sin(pkin(7));
t67 = sin(pkin(6));
t82 = t67 * t66;
t85 = (-t65 * t64 + t68 * t81) * t70 - t69 * t82;
t84 = t65 * t71;
t83 = t66 * t71;
t75 = cos(qJ(3));
t80 = t70 * t75;
t79 = t65 * t82;
t61 = -t69 * t64 - t68 * t84;
t62 = -t64 * t84 + t69 * t68;
t73 = sin(qJ(3));
t53 = t62 * t75 + (t61 * t70 + t79) * t73;
t57 = t65 * t67 * t70 - t61 * t66;
t72 = sin(qJ(4));
t74 = cos(qJ(4));
t44 = t53 * t74 + t57 * t72;
t42 = 0.1e1 / t44 ^ 2;
t43 = t53 * t72 - t57 * t74;
t77 = t43 ^ 2 * t42 + 0.1e1;
t60 = t64 * t81 + t65 * t68;
t56 = t73 * t83 + (t68 * t70 * t73 + t64 * t75) * t67;
t55 = -t75 * t83 + (t64 * t73 - t68 * t80) * t67;
t54 = 0.1e1 / t55 ^ 2;
t52 = -t61 * t80 + t62 * t73 - t75 * t79;
t51 = t60 * t75 + t85 * t73;
t49 = t60 * t73 - t85 * t75;
t48 = atan2(-t49, t55);
t46 = cos(t48);
t45 = sin(t48);
t41 = 0.1e1 / t77;
t40 = -t45 * t49 + t46 * t55;
t39 = 0.1e1 / t40 ^ 2;
t37 = (-t51 / t55 + t56 * t49 * t54) / (t49 ^ 2 * t54 + 0.1e1);
t1 = [0, 0, t37, 0, 0, 0; 0, 0 (t53 / t40 - (-t45 * t51 + t46 * t56 + (-t45 * t55 - t46 * t49) * t37) * t52 * t39) / (t52 ^ 2 * t39 + 0.1e1) 0, 0, 0; 0, 0 (-t72 / t44 + t74 * t43 * t42) * t52 * t41, t77 * t41, 0, 0;];
Ja_rot  = t1;
