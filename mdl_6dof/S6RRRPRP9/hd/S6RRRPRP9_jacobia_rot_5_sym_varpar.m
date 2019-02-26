% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP9_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:13:42
% EndTime: 2019-02-26 22:13:43
% DurationCPUTime: 0.12s
% Computational Cost: add. (141->25), mult. (425->68), div. (58->9), fcn. (611->11), ass. (0->40)
t60 = sin(qJ(1));
t62 = cos(qJ(3));
t63 = cos(qJ(2));
t58 = sin(qJ(3));
t64 = cos(qJ(1));
t67 = t64 * t58;
t47 = -t60 * t62 + t63 * t67;
t66 = t64 * t62;
t48 = t60 * t58 + t63 * t66;
t57 = sin(qJ(5));
t61 = cos(qJ(5));
t40 = t47 * t57 + t48 * t61;
t38 = 0.1e1 / t40 ^ 2;
t39 = -t47 * t61 + t48 * t57;
t75 = t38 * t39;
t74 = t39 ^ 2 * t38;
t59 = sin(qJ(2));
t69 = t60 * t59;
t52 = atan2(t69, t63);
t49 = sin(t52);
t50 = cos(t52);
t43 = t49 * t69 + t50 * t63;
t42 = 0.1e1 / t43 ^ 2;
t73 = t42 * t64 ^ 2;
t54 = t59 ^ 2;
t72 = t54 / t63 ^ 2;
t71 = t59 * t64;
t51 = 0.1e1 / (t60 ^ 2 * t72 + 0.1e1);
t70 = t60 * t51;
t68 = t60 * t63;
t65 = 0.1e1 + t74;
t55 = 0.1e1 / t63;
t46 = -t62 * t68 + t67;
t45 = -t58 * t68 - t66;
t44 = (0.1e1 + t72) * t70;
t41 = 0.1e1 / t43;
t37 = 0.1e1 / t40;
t36 = 0.1e1 / (t54 * t73 + 0.1e1);
t35 = 0.1e1 / t65;
t1 = [t55 * t51 * t71, t44, 0, 0, 0, 0; (t41 * t69 + (t50 * t54 * t55 * t70 + (-t51 + 0.1e1) * t59 * t49) * t59 * t73) * t36 (-t63 * t41 + (t49 * t68 - t50 * t59 + (-t49 * t63 + t50 * t69) * t44) * t59 * t42) * t64 * t36, 0, 0, 0, 0; ((-t45 * t61 + t46 * t57) * t37 - (t45 * t57 + t46 * t61) * t75) * t35 ((-t57 * t62 + t58 * t61) * t37 - (-t57 * t58 - t61 * t62) * t75) * t35 * t71 (-t37 * t40 - t74) * t35, 0, t65 * t35, 0;];
Ja_rot  = t1;
