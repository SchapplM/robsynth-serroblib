% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR7
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

function Ja_rot = S6RPPRRR7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:38:11
% EndTime: 2019-02-26 20:38:12
% DurationCPUTime: 0.11s
% Computational Cost: add. (519->20), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
t65 = pkin(10) + qJ(4) + qJ(5);
t63 = sin(t65);
t64 = cos(t65);
t69 = cos(qJ(1));
t73 = t69 * t64;
t58 = atan2(-t73, t63);
t52 = sin(t58);
t53 = cos(t58);
t49 = -t52 * t73 + t53 * t63;
t48 = 0.1e1 / t49 ^ 2;
t67 = sin(qJ(1));
t81 = t48 * t67 ^ 2;
t80 = t52 * t63;
t66 = sin(qJ(6));
t72 = t69 * t66;
t68 = cos(qJ(6));
t74 = t67 * t68;
t57 = t63 * t74 + t72;
t55 = 0.1e1 / t57 ^ 2;
t71 = t69 * t68;
t75 = t67 * t66;
t56 = t63 * t75 - t71;
t79 = t55 * t56;
t62 = t64 ^ 2;
t77 = 0.1e1 / t63 ^ 2 * t62;
t59 = 0.1e1 / (t69 ^ 2 * t77 + 0.1e1);
t78 = t59 * t69;
t76 = t64 * t67;
t70 = t56 ^ 2 * t55 + 0.1e1;
t60 = 0.1e1 / t63;
t54 = 0.1e1 / t57;
t51 = 0.1e1 / t70;
t50 = (0.1e1 + t77) * t78;
t47 = 0.1e1 / t49;
t46 = 0.1e1 / (t62 * t81 + 0.1e1);
t45 = (t54 * t66 - t68 * t79) * t51 * t76;
t44 = (t63 * t47 + (t69 * t80 + t53 * t64 + (-t53 * t73 - t80) * t50) * t64 * t48) * t67 * t46;
t1 = [t60 * t59 * t76, 0, 0, t50, t50, 0; (-t47 * t73 + (-t53 * t60 * t62 * t78 + (-t59 + 0.1e1) * t64 * t52) * t64 * t81) * t46, 0, 0, t44, t44, 0; ((t63 * t72 + t74) * t54 - (t63 * t71 - t75) * t79) * t51, 0, 0, t45, t45, t70 * t51;];
Ja_rot  = t1;
