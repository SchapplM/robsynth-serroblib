% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:47:10
% EndTime: 2019-02-26 21:47:10
% DurationCPUTime: 0.11s
% Computational Cost: add. (317->22), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->38)
t65 = qJ(2) + pkin(10);
t62 = cos(t65);
t61 = sin(t65);
t67 = sin(qJ(1));
t74 = t67 * t61;
t56 = atan2(-t74, -t62);
t54 = sin(t56);
t55 = cos(t56);
t47 = -t54 * t74 - t55 * t62;
t46 = 0.1e1 / t47 ^ 2;
t68 = cos(qJ(1));
t80 = t46 * t68 ^ 2;
t66 = qJ(4) + qJ(5);
t64 = cos(t66);
t70 = t68 * t64;
t63 = sin(t66);
t73 = t67 * t63;
t53 = t62 * t70 + t73;
t51 = 0.1e1 / t53 ^ 2;
t71 = t68 * t63;
t72 = t67 * t64;
t52 = t62 * t71 - t72;
t79 = t51 * t52;
t78 = t54 * t62;
t58 = t61 ^ 2;
t77 = t58 / t62 ^ 2;
t76 = t61 * t68;
t57 = 0.1e1 / (t67 ^ 2 * t77 + 0.1e1);
t75 = t67 * t57;
t69 = t52 ^ 2 * t51 + 0.1e1;
t59 = 0.1e1 / t62;
t50 = 0.1e1 / t53;
t49 = 0.1e1 / t69;
t48 = (0.1e1 + t77) * t75;
t45 = 0.1e1 / t47;
t44 = 0.1e1 / (t58 * t80 + 0.1e1);
t43 = t69 * t49;
t1 = [t59 * t57 * t76, t48, 0, 0, 0, 0; (-t45 * t74 - (-t55 * t58 * t59 * t75 + (t57 - 0.1e1) * t61 * t54) * t61 * t80) * t44 (t62 * t45 - (-t67 * t78 + t55 * t61 + (-t55 * t74 + t78) * t48) * t61 * t46) * t68 * t44, 0, 0, 0, 0; ((-t62 * t73 - t70) * t50 - (-t62 * t72 + t71) * t79) * t49 (-t50 * t63 + t64 * t79) * t49 * t76, 0, t43, t43, 0;];
Ja_rot  = t1;
