% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:11:36
% EndTime: 2019-02-26 22:11:36
% DurationCPUTime: 0.14s
% Computational Cost: add. (1004->28), mult. (689->68), div. (139->11), fcn. (1046->9), ass. (0->42)
t77 = sin(qJ(2));
t92 = t77 ^ 2;
t72 = qJ(3) + pkin(10) + qJ(5);
t70 = sin(t72);
t71 = cos(t72);
t80 = cos(qJ(1));
t82 = t80 * t71;
t78 = sin(qJ(1));
t79 = cos(qJ(2));
t84 = t78 * t79;
t61 = t70 * t84 + t82;
t86 = t77 * t70;
t57 = atan2(-t61, t86);
t54 = sin(t57);
t55 = cos(t57);
t52 = -t54 * t61 + t55 * t86;
t51 = 0.1e1 / t52 ^ 2;
t83 = t80 * t70;
t64 = -t78 * t71 + t79 * t83;
t91 = t51 * t64;
t89 = t55 * t61;
t88 = t64 ^ 2 * t51;
t68 = 0.1e1 / t70;
t74 = 0.1e1 / t77;
t87 = t68 * t74;
t85 = t77 * t80;
t65 = t78 * t70 + t79 * t82;
t60 = 0.1e1 / t65 ^ 2;
t81 = t80 ^ 2 * t92 * t60;
t75 = 0.1e1 / t92;
t69 = 0.1e1 / t70 ^ 2;
t63 = t71 * t84 - t83;
t59 = 0.1e1 / t65;
t58 = 0.1e1 / (0.1e1 + t81);
t56 = 0.1e1 / (t61 ^ 2 * t75 * t69 + 0.1e1);
t53 = t64 * t60 * t58 * t85;
t50 = 0.1e1 / t52;
t49 = (t61 * t68 * t75 * t79 + t78) * t56;
t48 = 0.1e1 / (0.1e1 + t88);
t47 = (t61 * t69 * t71 - t63 * t68) * t74 * t56;
t46 = (t65 * t50 - (t55 * t77 * t71 - t54 * t63 + (-t54 * t86 - t89) * t47) * t91) * t48;
t1 = [-t64 * t56 * t87, t49, t47, 0, t47, 0; (-t61 * t50 - (-t54 + (t87 * t89 + t54) * t56) * t88) * t48 (t49 * t89 * t91 + (-t50 * t85 - (t55 * t79 + (-t49 + t78) * t54 * t77) * t91) * t70) * t48, t46, 0, t46, 0; (-t60 * t63 * t80 + t59 * t78) * t77 * t58 (-t59 * t79 * t80 - t71 * t81) * t58, -t53, 0, -t53, 0;];
Ja_rot  = t1;
