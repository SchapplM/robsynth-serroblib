% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRR8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:38:39
% EndTime: 2019-02-26 20:38:39
% DurationCPUTime: 0.11s
% Computational Cost: add. (294->21), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->38)
t60 = pkin(10) + qJ(4);
t56 = sin(t60);
t57 = cos(t60);
t63 = cos(qJ(1));
t67 = t63 * t57;
t51 = atan2(-t67, t56);
t49 = sin(t51);
t50 = cos(t51);
t42 = -t49 * t67 + t50 * t56;
t41 = 0.1e1 / t42 ^ 2;
t62 = sin(qJ(1));
t75 = t41 * t62 ^ 2;
t61 = qJ(5) + qJ(6);
t58 = sin(t61);
t66 = t63 * t58;
t59 = cos(t61);
t69 = t62 * t59;
t48 = t56 * t69 + t66;
t46 = 0.1e1 / t48 ^ 2;
t65 = t63 * t59;
t70 = t62 * t58;
t47 = t56 * t70 - t65;
t74 = t46 * t47;
t73 = t49 * t56;
t55 = t57 ^ 2;
t72 = 0.1e1 / t56 ^ 2 * t55;
t71 = t57 * t62;
t52 = 0.1e1 / (t63 ^ 2 * t72 + 0.1e1);
t68 = t63 * t52;
t64 = t47 ^ 2 * t46 + 0.1e1;
t53 = 0.1e1 / t56;
t45 = 0.1e1 / t48;
t44 = 0.1e1 / t64;
t43 = (0.1e1 + t72) * t68;
t40 = 0.1e1 / t42;
t39 = 0.1e1 / (t55 * t75 + 0.1e1);
t38 = t64 * t44;
t1 = [t53 * t52 * t71, 0, 0, t43, 0, 0; (-t40 * t67 + (-t50 * t53 * t55 * t68 + (-t52 + 0.1e1) * t57 * t49) * t57 * t75) * t39, 0, 0 (t56 * t40 + (t63 * t73 + t50 * t57 + (-t50 * t67 - t73) * t43) * t57 * t41) * t62 * t39, 0, 0; ((t56 * t66 + t69) * t45 - (t56 * t65 - t70) * t74) * t44, 0, 0 (t45 * t58 - t59 * t74) * t44 * t71, t38, t38;];
Ja_rot  = t1;
