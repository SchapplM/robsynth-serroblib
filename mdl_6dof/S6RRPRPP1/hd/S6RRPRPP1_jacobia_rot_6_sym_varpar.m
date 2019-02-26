% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPP1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:51
% EndTime: 2019-02-26 21:34:52
% DurationCPUTime: 0.14s
% Computational Cost: add. (640->28), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->42)
t68 = qJ(2) + pkin(9);
t64 = sin(t68);
t84 = t64 ^ 2;
t66 = cos(t68);
t67 = qJ(4) + pkin(10);
t65 = cos(t67);
t71 = cos(qJ(1));
t73 = t71 * t65;
t63 = sin(t67);
t70 = sin(qJ(1));
t76 = t70 * t63;
t51 = t66 * t76 + t73;
t78 = t64 * t63;
t47 = atan2(-t51, t78);
t44 = sin(t47);
t45 = cos(t47);
t43 = -t44 * t51 + t45 * t78;
t42 = 0.1e1 / t43 ^ 2;
t74 = t71 * t63;
t75 = t70 * t65;
t54 = t66 * t74 - t75;
t83 = t42 * t54;
t81 = t45 * t51;
t80 = t54 ^ 2 * t42;
t58 = 0.1e1 / t63;
t61 = 0.1e1 / t64;
t79 = t58 * t61;
t77 = t64 * t71;
t55 = t66 * t73 + t76;
t50 = 0.1e1 / t55 ^ 2;
t72 = t71 ^ 2 * t84 * t50;
t62 = 0.1e1 / t84;
t59 = 0.1e1 / t63 ^ 2;
t53 = t66 * t75 - t74;
t49 = 0.1e1 / t55;
t48 = 0.1e1 / (0.1e1 + t72);
t46 = 0.1e1 / (t51 ^ 2 * t62 * t59 + 0.1e1);
t41 = 0.1e1 / t43;
t40 = (t51 * t58 * t62 * t66 + t70) * t46;
t39 = 0.1e1 / (0.1e1 + t80);
t38 = (t51 * t59 * t65 - t53 * t58) * t61 * t46;
t1 = [-t54 * t46 * t79, t40, 0, t38, 0, 0; (-t51 * t41 - (-t44 + (t79 * t81 + t44) * t46) * t80) * t39 (t40 * t81 * t83 + (-t41 * t77 - (t45 * t66 + (-t40 + t70) * t64 * t44) * t83) * t63) * t39, 0 (t55 * t41 - (t45 * t64 * t65 - t44 * t53 + (-t44 * t78 - t81) * t38) * t83) * t39, 0, 0; (-t50 * t53 * t71 + t49 * t70) * t64 * t48 (-t49 * t66 * t71 - t65 * t72) * t48, 0, -t54 * t50 * t48 * t77, 0, 0;];
Ja_rot  = t1;
