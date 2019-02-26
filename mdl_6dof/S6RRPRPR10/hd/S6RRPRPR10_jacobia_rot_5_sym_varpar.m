% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR10_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:00
% EndTime: 2019-02-26 21:43:00
% DurationCPUTime: 0.15s
% Computational Cost: add. (727->32), mult. (1022->76), div. (77->9), fcn. (1479->11), ass. (0->49)
t71 = cos(pkin(6));
t72 = sin(qJ(2));
t75 = cos(qJ(1));
t78 = t75 * t72;
t73 = sin(qJ(1));
t74 = cos(qJ(2));
t79 = t73 * t74;
t62 = t71 * t78 + t79;
t69 = pkin(11) + qJ(4);
t67 = sin(t69);
t68 = cos(t69);
t70 = sin(pkin(6));
t81 = t70 * t75;
t50 = t62 * t67 + t68 * t81;
t84 = t70 * t72;
t57 = t67 * t84 - t71 * t68;
t48 = atan2(-t50, t57);
t45 = sin(t48);
t46 = cos(t48);
t44 = -t45 * t50 + t46 * t57;
t43 = 0.1e1 / t44 ^ 2;
t77 = t75 * t74;
t80 = t73 * t72;
t64 = -t71 * t80 + t77;
t83 = t70 * t73;
t53 = t64 * t67 - t68 * t83;
t89 = t43 * t53;
t88 = t46 * t50;
t56 = 0.1e1 / t57 ^ 2;
t87 = t50 * t56;
t86 = t53 ^ 2 * t43;
t54 = t64 * t68 + t67 * t83;
t63 = t71 * t79 + t78;
t60 = 0.1e1 / t63 ^ 2;
t85 = t54 * t60;
t82 = t70 * t74;
t52 = t62 * t68 - t67 * t81;
t76 = -t45 * t57 - t88;
t61 = t71 * t77 - t80;
t59 = 0.1e1 / t63;
t58 = t71 * t67 + t68 * t84;
t55 = 0.1e1 / t57;
t49 = 0.1e1 / (t54 ^ 2 * t60 + 0.1e1);
t47 = 0.1e1 / (t50 ^ 2 * t56 + 0.1e1);
t42 = 0.1e1 / t44;
t41 = 0.1e1 / (0.1e1 + t86);
t40 = (-t55 * t61 + t82 * t87) * t67 * t47;
t39 = (-t52 * t55 + t58 * t87) * t47;
t1 = [-t53 * t55 * t47, t40, 0, t39, 0, 0; (-t50 * t42 - (-t45 + (t55 * t88 + t45) * t47) * t86) * t41 (-t63 * t67 * t42 - ((-t45 * t61 + t46 * t82) * t67 + t76 * t40) * t89) * t41, 0 (t54 * t42 - (t39 * t76 - t45 * t52 + t46 * t58) * t89) * t41, 0, 0; (-t52 * t59 - t61 * t85) * t49 (-t59 * t63 * t68 - t64 * t85) * t49, 0, -t53 * t59 * t49, 0, 0;];
Ja_rot  = t1;
