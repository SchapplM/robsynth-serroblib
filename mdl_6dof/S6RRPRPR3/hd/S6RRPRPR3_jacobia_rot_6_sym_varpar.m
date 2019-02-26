% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:38:51
% EndTime: 2019-02-26 21:38:51
% DurationCPUTime: 0.08s
% Computational Cost: add. (379->22), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->38)
t65 = qJ(2) + pkin(10);
t63 = cos(t65);
t62 = sin(t65);
t66 = sin(qJ(1));
t71 = t66 * t62;
t55 = atan2(-t71, -t63);
t53 = sin(t55);
t54 = cos(t55);
t47 = -t53 * t71 - t54 * t63;
t46 = 0.1e1 / t47 ^ 2;
t67 = cos(qJ(1));
t79 = t46 * t67 ^ 2;
t64 = qJ(4) + pkin(11) + qJ(6);
t58 = cos(t64);
t69 = t67 * t58;
t57 = sin(t64);
t73 = t66 * t57;
t52 = t63 * t69 + t73;
t50 = 0.1e1 / t52 ^ 2;
t70 = t67 * t57;
t72 = t66 * t58;
t51 = t63 * t70 - t72;
t78 = t50 * t51;
t77 = t53 * t63;
t59 = t62 ^ 2;
t76 = t59 / t63 ^ 2;
t75 = t62 * t67;
t56 = 0.1e1 / (t66 ^ 2 * t76 + 0.1e1);
t74 = t66 * t56;
t68 = t51 ^ 2 * t50 + 0.1e1;
t60 = 0.1e1 / t63;
t49 = 0.1e1 / t52;
t48 = (0.1e1 + t76) * t74;
t45 = 0.1e1 / t47;
t44 = 0.1e1 / t68;
t43 = 0.1e1 / (t59 * t79 + 0.1e1);
t42 = t68 * t44;
t1 = [t60 * t56 * t75, t48, 0, 0, 0, 0; (-t45 * t71 - (-t54 * t59 * t60 * t74 + (t56 - 0.1e1) * t62 * t53) * t62 * t79) * t43 (t63 * t45 - (-t66 * t77 + t54 * t62 + (-t54 * t71 + t77) * t48) * t62 * t46) * t67 * t43, 0, 0, 0, 0; ((-t63 * t73 - t69) * t49 - (-t63 * t72 + t70) * t78) * t44 (-t49 * t57 + t58 * t78) * t44 * t75, 0, t42, 0, t42;];
Ja_rot  = t1;
