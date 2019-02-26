% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:15:25
% EndTime: 2019-02-26 21:15:25
% DurationCPUTime: 0.08s
% Computational Cost: add. (517->22), mult. (332->57), div. (79->9), fcn. (489->9), ass. (0->38)
t67 = qJ(1) + pkin(11);
t64 = cos(t67);
t81 = t64 ^ 2;
t68 = qJ(3) + qJ(4);
t66 = cos(t68);
t63 = sin(t67);
t65 = sin(t68);
t75 = t63 * t65;
t58 = atan2(-t75, -t66);
t56 = sin(t58);
t57 = cos(t58);
t49 = -t56 * t75 - t57 * t66;
t48 = 0.1e1 / t49 ^ 2;
t80 = t48 * t65;
t69 = sin(qJ(5));
t70 = cos(qJ(5));
t72 = t66 * t70;
t55 = t63 * t69 + t64 * t72;
t53 = 0.1e1 / t55 ^ 2;
t73 = t66 * t69;
t54 = -t63 * t70 + t64 * t73;
t79 = t53 * t54;
t78 = t56 * t66;
t60 = t65 ^ 2;
t77 = t60 / t66 ^ 2;
t59 = 0.1e1 / (t63 ^ 2 * t77 + 0.1e1);
t76 = t63 * t59;
t74 = t64 * t65;
t71 = t54 ^ 2 * t53 + 0.1e1;
t61 = 0.1e1 / t66;
t52 = 0.1e1 / t55;
t51 = 0.1e1 / t71;
t50 = (0.1e1 + t77) * t76;
t47 = 0.1e1 / t49;
t46 = 0.1e1 / (t81 * t60 * t48 + 0.1e1);
t45 = (-t52 * t69 + t70 * t79) * t51 * t74;
t44 = (t66 * t47 - (-t63 * t78 + t57 * t65 + (-t57 * t75 + t78) * t50) * t80) * t64 * t46;
t1 = [t61 * t59 * t74, 0, t50, t50, 0, 0; (-t47 * t75 - (-t57 * t60 * t61 * t76 + (t59 - 0.1e1) * t65 * t56) * t81 * t80) * t46, 0, t44, t44, 0, 0; ((-t63 * t73 - t64 * t70) * t52 - (-t63 * t72 + t64 * t69) * t79) * t51, 0, t45, t45, t71 * t51, 0;];
Ja_rot  = t1;
