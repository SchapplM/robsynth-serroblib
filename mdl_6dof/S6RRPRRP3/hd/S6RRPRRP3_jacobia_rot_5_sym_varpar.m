% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_rot = S6RRPRRP3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:47:14
% EndTime: 2019-02-26 21:47:14
% DurationCPUTime: 0.07s
% Computational Cost: add. (317->22), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->38)
t64 = qJ(2) + pkin(10);
t61 = cos(t64);
t60 = sin(t64);
t66 = sin(qJ(1));
t73 = t66 * t60;
t55 = atan2(-t73, -t61);
t53 = sin(t55);
t54 = cos(t55);
t46 = -t53 * t73 - t54 * t61;
t45 = 0.1e1 / t46 ^ 2;
t67 = cos(qJ(1));
t79 = t45 * t67 ^ 2;
t65 = qJ(4) + qJ(5);
t63 = cos(t65);
t69 = t67 * t63;
t62 = sin(t65);
t72 = t66 * t62;
t52 = t61 * t69 + t72;
t50 = 0.1e1 / t52 ^ 2;
t70 = t67 * t62;
t71 = t66 * t63;
t51 = t61 * t70 - t71;
t78 = t50 * t51;
t77 = t53 * t61;
t57 = t60 ^ 2;
t76 = t57 / t61 ^ 2;
t75 = t60 * t67;
t56 = 0.1e1 / (t66 ^ 2 * t76 + 0.1e1);
t74 = t66 * t56;
t68 = t51 ^ 2 * t50 + 0.1e1;
t58 = 0.1e1 / t61;
t49 = 0.1e1 / t52;
t48 = 0.1e1 / t68;
t47 = (0.1e1 + t76) * t74;
t44 = 0.1e1 / t46;
t43 = 0.1e1 / (t57 * t79 + 0.1e1);
t42 = t68 * t48;
t1 = [t58 * t56 * t75, t47, 0, 0, 0, 0; (-t44 * t73 - (-t54 * t57 * t58 * t74 + (t56 - 0.1e1) * t60 * t53) * t60 * t79) * t43 (t61 * t44 - (-t66 * t77 + t54 * t60 + (-t54 * t73 + t77) * t47) * t60 * t45) * t67 * t43, 0, 0, 0, 0; ((-t61 * t72 - t69) * t49 - (-t61 * t71 + t70) * t78) * t48 (-t49 * t62 + t63 * t78) * t48 * t75, 0, t42, t42, 0;];
Ja_rot  = t1;
