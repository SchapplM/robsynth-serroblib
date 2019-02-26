% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR3
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
% Datum: 2019-02-26 21:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:16:00
% EndTime: 2019-02-26 21:16:00
% DurationCPUTime: 0.11s
% Computational Cost: add. (418->22), mult. (278->57), div. (62->9), fcn. (402->9), ass. (0->36)
t60 = qJ(1) + pkin(11);
t58 = cos(t60);
t75 = t58 ^ 2;
t65 = cos(qJ(3));
t57 = sin(t60);
t64 = sin(qJ(3));
t71 = t57 * t64;
t53 = atan2(-t71, -t65);
t51 = sin(t53);
t52 = cos(t53);
t45 = -t51 * t71 - t52 * t65;
t44 = 0.1e1 / t45 ^ 2;
t74 = t44 * t64;
t59 = qJ(4) + qJ(5) + qJ(6);
t55 = sin(t59);
t56 = cos(t59);
t68 = t58 * t65;
t50 = t57 * t55 + t56 * t68;
t48 = 0.1e1 / t50 ^ 2;
t49 = t55 * t68 - t57 * t56;
t73 = t48 * t49;
t61 = t64 ^ 2;
t67 = t61 / t65 ^ 2;
t54 = 0.1e1 / (t57 ^ 2 * t67 + 0.1e1);
t72 = t57 * t54;
t70 = t57 * t65;
t69 = t58 * t64;
t66 = t49 ^ 2 * t48 + 0.1e1;
t62 = 0.1e1 / t65;
t47 = 0.1e1 / t50;
t46 = (0.1e1 + t67) * t72;
t43 = 0.1e1 / t45;
t42 = 0.1e1 / (t75 * t61 * t44 + 0.1e1);
t41 = 0.1e1 / t66;
t40 = t66 * t41;
t1 = [t62 * t54 * t69, 0, t46, 0, 0, 0; (-t43 * t71 - (-t52 * t61 * t62 * t72 + (t54 - 0.1e1) * t64 * t51) * t75 * t74) * t42, 0 (t65 * t43 - (-t51 * t70 + t52 * t64 + (t51 * t65 - t52 * t71) * t46) * t74) * t58 * t42, 0, 0, 0; ((-t55 * t70 - t58 * t56) * t47 - (t58 * t55 - t56 * t70) * t73) * t41, 0 (-t47 * t55 + t56 * t73) * t41 * t69, t40, t40, t40;];
Ja_rot  = t1;
