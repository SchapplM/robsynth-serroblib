% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR6
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
% Datum: 2019-02-26 21:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:17:36
% EndTime: 2019-02-26 21:17:36
% DurationCPUTime: 0.11s
% Computational Cost: add. (424->22), mult. (278->54), div. (62->9), fcn. (402->9), ass. (0->38)
t68 = pkin(11) + qJ(3);
t66 = cos(t68);
t65 = sin(t68);
t69 = sin(qJ(1));
t74 = t69 * t65;
t58 = atan2(-t74, -t66);
t56 = sin(t58);
t57 = cos(t58);
t50 = -t56 * t74 - t57 * t66;
t48 = 0.1e1 / t50 ^ 2;
t70 = cos(qJ(1));
t82 = t48 * t70 ^ 2;
t67 = qJ(4) + qJ(5) + qJ(6);
t61 = cos(t67);
t72 = t70 * t61;
t60 = sin(t67);
t76 = t69 * t60;
t55 = t66 * t72 + t76;
t53 = 0.1e1 / t55 ^ 2;
t73 = t70 * t60;
t75 = t69 * t61;
t54 = t66 * t73 - t75;
t81 = t53 * t54;
t80 = t56 * t66;
t62 = t65 ^ 2;
t79 = t62 / t66 ^ 2;
t78 = t65 * t70;
t59 = 0.1e1 / (t69 ^ 2 * t79 + 0.1e1);
t77 = t69 * t59;
t71 = t54 ^ 2 * t53 + 0.1e1;
t63 = 0.1e1 / t66;
t52 = 0.1e1 / t55;
t51 = (0.1e1 + t79) * t77;
t49 = 0.1e1 / t71;
t47 = 0.1e1 / t50;
t46 = 0.1e1 / (t62 * t82 + 0.1e1);
t45 = t71 * t49;
t1 = [t63 * t59 * t78, 0, t51, 0, 0, 0; (-t47 * t74 - (-t57 * t62 * t63 * t77 + (t59 - 0.1e1) * t65 * t56) * t65 * t82) * t46, 0 (t66 * t47 - (-t69 * t80 + t57 * t65 + (-t57 * t74 + t80) * t51) * t65 * t48) * t70 * t46, 0, 0, 0; ((-t66 * t76 - t72) * t52 - (-t66 * t75 + t73) * t81) * t49, 0 (-t52 * t60 + t61 * t81) * t49 * t78, t45, t45, t45;];
Ja_rot  = t1;
