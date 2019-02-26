% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP2
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
% Datum: 2019-02-26 21:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:08:37
% EndTime: 2019-02-26 21:08:37
% DurationCPUTime: 0.10s
% Computational Cost: add. (304->22), mult. (251->57), div. (57->9), fcn. (367->9), ass. (0->37)
t55 = qJ(1) + pkin(10);
t52 = cos(t55);
t72 = t52 ^ 2;
t61 = cos(qJ(3));
t51 = sin(t55);
t60 = sin(qJ(3));
t67 = t51 * t60;
t49 = atan2(-t67, -t61);
t47 = sin(t49);
t48 = cos(t49);
t41 = -t47 * t67 - t48 * t61;
t40 = 0.1e1 / t41 ^ 2;
t71 = t40 * t60;
t59 = qJ(4) + qJ(5);
t53 = sin(t59);
t54 = cos(t59);
t64 = t54 * t61;
t46 = t51 * t53 + t52 * t64;
t44 = 0.1e1 / t46 ^ 2;
t65 = t53 * t61;
t45 = -t51 * t54 + t52 * t65;
t70 = t44 * t45;
t69 = t47 * t61;
t56 = t60 ^ 2;
t63 = t56 / t61 ^ 2;
t50 = 0.1e1 / (t51 ^ 2 * t63 + 0.1e1);
t68 = t51 * t50;
t66 = t52 * t60;
t62 = t45 ^ 2 * t44 + 0.1e1;
t57 = 0.1e1 / t61;
t43 = 0.1e1 / t46;
t42 = (0.1e1 + t63) * t68;
t39 = 0.1e1 / t41;
t38 = 0.1e1 / t62;
t37 = 0.1e1 / (t72 * t56 * t40 + 0.1e1);
t36 = t62 * t38;
t1 = [t57 * t50 * t66, 0, t42, 0, 0, 0; (-t39 * t67 - (-t48 * t56 * t57 * t68 + (t50 - 0.1e1) * t60 * t47) * t72 * t71) * t37, 0 (t61 * t39 - (-t51 * t69 + t48 * t60 + (-t48 * t67 + t69) * t42) * t71) * t52 * t37, 0, 0, 0; ((-t51 * t65 - t52 * t54) * t43 - (-t51 * t64 + t52 * t53) * t70) * t38, 0 (-t43 * t53 + t54 * t70) * t38 * t66, t36, t36, 0;];
Ja_rot  = t1;
