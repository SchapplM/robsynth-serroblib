% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:01
% EndTime: 2019-02-26 21:02:02
% DurationCPUTime: 0.12s
% Computational Cost: add. (334->26), mult. (425->70), div. (58->9), fcn. (611->11), ass. (0->42)
t56 = qJ(1) + pkin(10);
t55 = cos(t56);
t77 = t55 ^ 2;
t54 = sin(t56);
t64 = cos(qJ(4));
t61 = sin(qJ(4));
t65 = cos(qJ(3));
t68 = t61 * t65;
t47 = -t54 * t64 + t55 * t68;
t67 = t64 * t65;
t48 = t54 * t61 + t55 * t67;
t60 = sin(qJ(6));
t63 = cos(qJ(6));
t40 = t47 * t60 + t48 * t63;
t38 = 0.1e1 / t40 ^ 2;
t39 = -t47 * t63 + t48 * t60;
t76 = t38 * t39;
t75 = t39 ^ 2 * t38;
t62 = sin(qJ(3));
t71 = t54 * t62;
t52 = atan2(t71, t65);
t49 = sin(t52);
t50 = cos(t52);
t44 = t49 * t71 + t50 * t65;
t43 = 0.1e1 / t44 ^ 2;
t74 = t43 * t62;
t73 = t49 * t65;
t57 = t62 ^ 2;
t69 = t57 / t65 ^ 2;
t51 = 0.1e1 / (t54 ^ 2 * t69 + 0.1e1);
t72 = t54 * t51;
t70 = t55 * t62;
t66 = 0.1e1 + t75;
t58 = 0.1e1 / t65;
t46 = -t54 * t67 + t55 * t61;
t45 = -t54 * t68 - t55 * t64;
t42 = 0.1e1 / t44;
t41 = (0.1e1 + t69) * t72;
t37 = 0.1e1 / t40;
t36 = 0.1e1 / (t77 * t57 * t43 + 0.1e1);
t35 = 0.1e1 / t66;
t1 = [t58 * t51 * t70, 0, t41, 0, 0, 0; (t42 * t71 + (t50 * t57 * t58 * t72 + (-t51 + 0.1e1) * t62 * t49) * t77 * t74) * t36, 0 (-t65 * t42 + (t54 * t73 - t50 * t62 + (t50 * t71 - t73) * t41) * t74) * t55 * t36, 0, 0, 0; ((-t45 * t63 + t46 * t60) * t37 - (t45 * t60 + t46 * t63) * t76) * t35, 0 ((-t60 * t64 + t61 * t63) * t37 - (-t60 * t61 - t63 * t64) * t76) * t35 * t70 (-t40 * t37 - t75) * t35, 0, t66 * t35;];
Ja_rot  = t1;
