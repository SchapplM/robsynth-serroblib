% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR9_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:14
% EndTime: 2019-02-26 21:19:14
% DurationCPUTime: 0.11s
% Computational Cost: add. (258->20), mult. (278->55), div. (62->9), fcn. (402->9), ass. (0->38)
t61 = sin(qJ(1));
t76 = t61 ^ 2;
t60 = sin(qJ(3));
t62 = cos(qJ(3));
t63 = cos(qJ(1));
t65 = t63 * t62;
t52 = atan2(-t65, t60);
t50 = sin(t52);
t51 = cos(t52);
t44 = -t50 * t65 + t51 * t60;
t43 = 0.1e1 / t44 ^ 2;
t75 = t43 * t62;
t56 = qJ(4) + qJ(5) + qJ(6);
t54 = sin(t56);
t67 = t63 * t54;
t55 = cos(t56);
t70 = t61 * t55;
t49 = t60 * t70 + t67;
t47 = 0.1e1 / t49 ^ 2;
t66 = t63 * t55;
t71 = t61 * t54;
t48 = t60 * t71 - t66;
t74 = t47 * t48;
t73 = t50 * t60;
t59 = t62 ^ 2;
t72 = 0.1e1 / t60 ^ 2 * t59;
t69 = t61 * t62;
t53 = 0.1e1 / (t63 ^ 2 * t72 + 0.1e1);
t68 = t63 * t53;
t64 = t48 ^ 2 * t47 + 0.1e1;
t57 = 0.1e1 / t60;
t46 = 0.1e1 / t49;
t45 = (0.1e1 + t72) * t68;
t42 = 0.1e1 / t44;
t41 = 0.1e1 / (t76 * t59 * t43 + 0.1e1);
t40 = 0.1e1 / t64;
t39 = t64 * t40;
t1 = [t57 * t53 * t69, 0, t45, 0, 0, 0; (-t42 * t65 + (-t51 * t57 * t59 * t68 + (-t53 + 0.1e1) * t62 * t50) * t76 * t75) * t41, 0 (t60 * t42 + (t63 * t73 + t51 * t62 + (-t51 * t65 - t73) * t45) * t75) * t61 * t41, 0, 0, 0; ((t60 * t67 + t70) * t46 - (t60 * t66 - t71) * t74) * t40, 0 (t46 * t54 - t55 * t74) * t40 * t69, t39, t39, t39;];
Ja_rot  = t1;
