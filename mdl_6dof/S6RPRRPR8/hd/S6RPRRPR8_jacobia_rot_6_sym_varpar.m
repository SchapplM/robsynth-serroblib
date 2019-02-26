% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:01
% EndTime: 2019-02-26 21:05:02
% DurationCPUTime: 0.10s
% Computational Cost: add. (220->20), mult. (251->55), div. (57->9), fcn. (367->9), ass. (0->38)
t57 = sin(qJ(1));
t72 = t57 ^ 2;
t56 = sin(qJ(3));
t58 = cos(qJ(3));
t59 = cos(qJ(1));
t61 = t59 * t58;
t48 = atan2(-t61, t56);
t46 = sin(t48);
t47 = cos(t48);
t40 = -t46 * t61 + t47 * t56;
t39 = 0.1e1 / t40 ^ 2;
t71 = t39 * t58;
t52 = qJ(4) + pkin(10) + qJ(6);
t50 = sin(t52);
t63 = t59 * t50;
t51 = cos(t52);
t66 = t57 * t51;
t45 = t56 * t66 + t63;
t43 = 0.1e1 / t45 ^ 2;
t62 = t59 * t51;
t67 = t57 * t50;
t44 = t56 * t67 - t62;
t70 = t43 * t44;
t69 = t46 * t56;
t55 = t58 ^ 2;
t68 = 0.1e1 / t56 ^ 2 * t55;
t65 = t57 * t58;
t49 = 0.1e1 / (t59 ^ 2 * t68 + 0.1e1);
t64 = t59 * t49;
t60 = t44 ^ 2 * t43 + 0.1e1;
t53 = 0.1e1 / t56;
t42 = 0.1e1 / t45;
t41 = (0.1e1 + t68) * t64;
t38 = 0.1e1 / t40;
t37 = 0.1e1 / (t72 * t55 * t39 + 0.1e1);
t36 = 0.1e1 / t60;
t35 = t60 * t36;
t1 = [t53 * t49 * t65, 0, t41, 0, 0, 0; (-t38 * t61 + (-t47 * t53 * t55 * t64 + (-t49 + 0.1e1) * t58 * t46) * t72 * t71) * t37, 0 (t56 * t38 + (t59 * t69 + t47 * t58 + (-t47 * t61 - t69) * t41) * t71) * t57 * t37, 0, 0, 0; ((t56 * t63 + t66) * t42 - (t56 * t62 - t67) * t70) * t36, 0 (t42 * t50 - t51 * t70) * t36 * t65, t35, 0, t35;];
Ja_rot  = t1;
