% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR6
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
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:47
% EndTime: 2019-02-26 22:18:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (357->21), mult. (278->55), div. (62->9), fcn. (402->9), ass. (0->35)
t65 = cos(qJ(2));
t63 = sin(qJ(2));
t64 = sin(qJ(1));
t71 = t64 * t63;
t55 = atan2(-t71, -t65);
t53 = sin(t55);
t54 = cos(t55);
t47 = -t53 * t71 - t54 * t65;
t46 = 0.1e1 / t47 ^ 2;
t66 = cos(qJ(1));
t76 = t46 * t66 ^ 2;
t59 = qJ(3) + pkin(11) + qJ(5) + qJ(6);
t57 = sin(t59);
t58 = cos(t59);
t68 = t66 * t58;
t52 = t64 * t57 + t65 * t68;
t50 = 0.1e1 / t52 ^ 2;
t69 = t66 * t57;
t51 = -t64 * t58 + t65 * t69;
t75 = t50 * t51;
t60 = t63 ^ 2;
t74 = t60 / t65 ^ 2;
t73 = t63 * t66;
t56 = 0.1e1 / (t64 ^ 2 * t74 + 0.1e1);
t72 = t64 * t56;
t70 = t64 * t65;
t67 = t51 ^ 2 * t50 + 0.1e1;
t61 = 0.1e1 / t65;
t49 = 0.1e1 / t52;
t48 = (0.1e1 + t74) * t72;
t45 = 0.1e1 / t47;
t44 = 0.1e1 / (t60 * t76 + 0.1e1);
t43 = 0.1e1 / t67;
t42 = t67 * t43;
t1 = [t61 * t56 * t73, t48, 0, 0, 0, 0; (-t45 * t71 - (-t54 * t60 * t61 * t72 + (t56 - 0.1e1) * t63 * t53) * t63 * t76) * t44 (t65 * t45 - (-t53 * t70 + t54 * t63 + (t53 * t65 - t54 * t71) * t48) * t63 * t46) * t66 * t44, 0, 0, 0, 0; ((-t57 * t70 - t68) * t49 - (-t58 * t70 + t69) * t75) * t43 (-t49 * t57 + t58 * t75) * t43 * t73, t42, 0, t42, t42;];
Ja_rot  = t1;
