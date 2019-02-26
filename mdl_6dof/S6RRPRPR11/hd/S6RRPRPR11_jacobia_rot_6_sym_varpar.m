% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR11_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:36
% EndTime: 2019-02-26 21:43:36
% DurationCPUTime: 0.10s
% Computational Cost: add. (221->21), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->37)
t58 = sin(qJ(2));
t59 = sin(qJ(1));
t60 = cos(qJ(2));
t66 = t59 * t60;
t50 = atan2(-t66, t58);
t48 = sin(t50);
t49 = cos(t50);
t42 = -t48 * t66 + t49 * t58;
t41 = 0.1e1 / t42 ^ 2;
t61 = cos(qJ(1));
t73 = t41 * t61 ^ 2;
t54 = qJ(4) + pkin(10) + qJ(6);
t52 = sin(t54);
t64 = t61 * t52;
t53 = cos(t54);
t67 = t59 * t53;
t47 = t58 * t64 + t67;
t45 = 0.1e1 / t47 ^ 2;
t63 = t61 * t53;
t68 = t59 * t52;
t46 = -t58 * t63 + t68;
t72 = t45 * t46;
t71 = t48 * t58;
t57 = t60 ^ 2;
t70 = 0.1e1 / t58 ^ 2 * t57;
t51 = 0.1e1 / (t59 ^ 2 * t70 + 0.1e1);
t69 = t59 * t51;
t65 = t60 * t61;
t62 = t46 ^ 2 * t45 + 0.1e1;
t55 = 0.1e1 / t58;
t44 = 0.1e1 / t47;
t43 = (0.1e1 + t70) * t69;
t40 = 0.1e1 / t42;
t39 = 0.1e1 / (t57 * t73 + 0.1e1);
t38 = 0.1e1 / t62;
t37 = t62 * t38;
t1 = [-t55 * t51 * t65, t43, 0, 0, 0, 0; (-t40 * t66 - (t49 * t55 * t57 * t69 + (t51 - 0.1e1) * t60 * t48) * t60 * t73) * t39 (-t58 * t40 - (t59 * t71 + t49 * t60 + (-t49 * t66 - t71) * t43) * t60 * t41) * t61 * t39, 0, 0, 0, 0; ((t58 * t67 + t64) * t44 - (-t58 * t68 + t63) * t72) * t38 (-t44 * t53 - t52 * t72) * t38 * t65, 0, t37, 0, t37;];
Ja_rot  = t1;
