% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR4_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:17:42
% EndTime: 2019-02-26 22:17:42
% DurationCPUTime: 0.11s
% Computational Cost: add. (341->21), mult. (305->53), div. (74->9), fcn. (454->9), ass. (0->37)
t62 = qJ(2) + qJ(3);
t61 = cos(t62);
t60 = sin(t62);
t65 = sin(qJ(1));
t71 = t65 * t60;
t55 = atan2(-t71, -t61);
t53 = sin(t55);
t54 = cos(t55);
t46 = -t53 * t71 - t54 * t61;
t45 = 0.1e1 / t46 ^ 2;
t66 = cos(qJ(1));
t77 = t45 * t66 ^ 2;
t64 = cos(pkin(11));
t67 = t66 * t64;
t63 = sin(pkin(11));
t70 = t65 * t63;
t52 = t61 * t67 + t70;
t50 = 0.1e1 / t52 ^ 2;
t68 = t66 * t63;
t69 = t65 * t64;
t51 = t61 * t68 - t69;
t76 = t50 * t51;
t75 = t53 * t61;
t57 = t60 ^ 2;
t74 = t57 / t61 ^ 2;
t73 = t60 * t66;
t56 = 0.1e1 / (t65 ^ 2 * t74 + 0.1e1);
t72 = t65 * t56;
t58 = 0.1e1 / t61;
t49 = 0.1e1 / t52;
t48 = 0.1e1 / (t51 ^ 2 * t50 + 0.1e1);
t47 = (0.1e1 + t74) * t72;
t44 = 0.1e1 / t46;
t43 = 0.1e1 / (t57 * t77 + 0.1e1);
t42 = (-t49 * t63 + t64 * t76) * t48 * t73;
t41 = (t61 * t44 - (-t65 * t75 + t54 * t60 + (-t54 * t71 + t75) * t47) * t60 * t45) * t66 * t43;
t1 = [t58 * t56 * t73, t47, t47, 0, 0, 0; (-t44 * t71 - (-t54 * t57 * t58 * t72 + (t56 - 0.1e1) * t60 * t53) * t60 * t77) * t43, t41, t41, 0, 0, 0; ((-t61 * t70 - t67) * t49 - (-t61 * t69 + t68) * t76) * t48, t42, t42, 0, 0, 0;];
Ja_rot  = t1;
