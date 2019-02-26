% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:10:50
% EndTime: 2019-02-26 21:10:50
% DurationCPUTime: 0.08s
% Computational Cost: add. (317->22), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->38)
t62 = pkin(10) + qJ(3);
t59 = cos(t62);
t58 = sin(t62);
t64 = sin(qJ(1));
t71 = t64 * t58;
t53 = atan2(-t71, -t59);
t51 = sin(t53);
t52 = cos(t53);
t44 = -t51 * t71 - t52 * t59;
t43 = 0.1e1 / t44 ^ 2;
t65 = cos(qJ(1));
t77 = t43 * t65 ^ 2;
t63 = qJ(4) + qJ(5);
t61 = cos(t63);
t67 = t65 * t61;
t60 = sin(t63);
t70 = t64 * t60;
t50 = t59 * t67 + t70;
t48 = 0.1e1 / t50 ^ 2;
t68 = t65 * t60;
t69 = t64 * t61;
t49 = t59 * t68 - t69;
t76 = t48 * t49;
t75 = t51 * t59;
t55 = t58 ^ 2;
t74 = t55 / t59 ^ 2;
t73 = t58 * t65;
t54 = 0.1e1 / (t64 ^ 2 * t74 + 0.1e1);
t72 = t64 * t54;
t66 = t49 ^ 2 * t48 + 0.1e1;
t56 = 0.1e1 / t59;
t47 = 0.1e1 / t50;
t46 = 0.1e1 / t66;
t45 = (0.1e1 + t74) * t72;
t42 = 0.1e1 / t44;
t41 = 0.1e1 / (t55 * t77 + 0.1e1);
t40 = t66 * t46;
t1 = [t56 * t54 * t73, 0, t45, 0, 0, 0; (-t42 * t71 - (-t52 * t55 * t56 * t72 + (t54 - 0.1e1) * t58 * t51) * t58 * t77) * t41, 0 (t59 * t42 - (-t64 * t75 + t52 * t58 + (-t52 * t71 + t75) * t45) * t58 * t43) * t65 * t41, 0, 0, 0; ((-t59 * t70 - t67) * t47 - (-t59 * t69 + t68) * t76) * t46, 0 (-t47 * t60 + t61 * t76) * t46 * t73, t40, t40, 0;];
Ja_rot  = t1;
