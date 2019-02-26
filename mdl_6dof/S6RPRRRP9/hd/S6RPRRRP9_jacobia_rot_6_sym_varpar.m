% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP9_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:12:30
% EndTime: 2019-02-26 21:12:30
% DurationCPUTime: 0.10s
% Computational Cost: add. (158->20), mult. (251->55), div. (57->9), fcn. (367->9), ass. (0->38)
t56 = sin(qJ(1));
t71 = t56 ^ 2;
t55 = sin(qJ(3));
t57 = cos(qJ(3));
t58 = cos(qJ(1));
t60 = t58 * t57;
t47 = atan2(-t60, t55);
t45 = sin(t47);
t46 = cos(t47);
t39 = -t45 * t60 + t46 * t55;
t38 = 0.1e1 / t39 ^ 2;
t70 = t38 * t57;
t54 = qJ(4) + qJ(5);
t49 = sin(t54);
t62 = t58 * t49;
t50 = cos(t54);
t65 = t56 * t50;
t44 = t55 * t65 + t62;
t42 = 0.1e1 / t44 ^ 2;
t61 = t58 * t50;
t66 = t56 * t49;
t43 = t55 * t66 - t61;
t69 = t42 * t43;
t68 = t45 * t55;
t53 = t57 ^ 2;
t67 = 0.1e1 / t55 ^ 2 * t53;
t64 = t56 * t57;
t48 = 0.1e1 / (t58 ^ 2 * t67 + 0.1e1);
t63 = t58 * t48;
t59 = t43 ^ 2 * t42 + 0.1e1;
t51 = 0.1e1 / t55;
t41 = 0.1e1 / t44;
t40 = (0.1e1 + t67) * t63;
t37 = 0.1e1 / t39;
t36 = 0.1e1 / t59;
t35 = 0.1e1 / (t71 * t53 * t38 + 0.1e1);
t34 = t59 * t36;
t1 = [t51 * t48 * t64, 0, t40, 0, 0, 0; (-t37 * t60 + (-t46 * t51 * t53 * t63 + (-t48 + 0.1e1) * t57 * t45) * t71 * t70) * t35, 0 (t55 * t37 + (t58 * t68 + t46 * t57 + (-t46 * t60 - t68) * t40) * t70) * t56 * t35, 0, 0, 0; ((t55 * t62 + t65) * t41 - (t55 * t61 - t66) * t69) * t36, 0 (t41 * t49 - t50 * t69) * t36 * t64, t34, t34, 0;];
Ja_rot  = t1;
