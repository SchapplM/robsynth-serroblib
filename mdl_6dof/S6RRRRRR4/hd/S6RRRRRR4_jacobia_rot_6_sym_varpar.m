% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:48:47
% EndTime: 2019-02-26 22:48:47
% DurationCPUTime: 0.08s
% Computational Cost: add. (409->21), mult. (305->55), div. (67->9), fcn. (437->9), ass. (0->35)
t68 = cos(qJ(2));
t66 = sin(qJ(2));
t67 = sin(qJ(1));
t74 = t67 * t66;
t58 = atan2(-t74, -t68);
t56 = sin(t58);
t57 = cos(t58);
t50 = -t56 * t74 - t57 * t68;
t49 = 0.1e1 / t50 ^ 2;
t69 = cos(qJ(1));
t79 = t49 * t69 ^ 2;
t62 = qJ(3) + qJ(4) + qJ(5) + qJ(6);
t60 = sin(t62);
t61 = cos(t62);
t71 = t69 * t61;
t55 = t67 * t60 + t68 * t71;
t53 = 0.1e1 / t55 ^ 2;
t72 = t69 * t60;
t54 = -t67 * t61 + t68 * t72;
t78 = t53 * t54;
t63 = t66 ^ 2;
t77 = t63 / t68 ^ 2;
t76 = t66 * t69;
t59 = 0.1e1 / (t67 ^ 2 * t77 + 0.1e1);
t75 = t67 * t59;
t73 = t67 * t68;
t70 = t54 ^ 2 * t53 + 0.1e1;
t64 = 0.1e1 / t68;
t52 = 0.1e1 / t55;
t51 = (0.1e1 + t77) * t75;
t48 = 0.1e1 / t50;
t47 = 0.1e1 / (t63 * t79 + 0.1e1);
t46 = 0.1e1 / t70;
t45 = t70 * t46;
t1 = [t64 * t59 * t76, t51, 0, 0, 0, 0; (-t48 * t74 - (-t57 * t63 * t64 * t75 + (t59 - 0.1e1) * t66 * t56) * t66 * t79) * t47 (t68 * t48 - (-t56 * t73 + t57 * t66 + (t56 * t68 - t57 * t74) * t51) * t66 * t49) * t69 * t47, 0, 0, 0, 0; ((-t60 * t73 - t71) * t52 - (-t61 * t73 + t72) * t78) * t46 (-t52 * t60 + t61 * t78) * t46 * t76, t45, t45, t45, t45;];
Ja_rot  = t1;
