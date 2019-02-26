% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:10:35
% EndTime: 2019-02-26 22:10:35
% DurationCPUTime: 0.11s
% Computational Cost: add. (422->22), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->39)
t66 = qJ(2) + qJ(3);
t64 = cos(t66);
t63 = sin(t66);
t67 = sin(qJ(1));
t72 = t67 * t63;
t56 = atan2(-t72, -t64);
t54 = sin(t56);
t55 = cos(t56);
t47 = -t54 * t72 - t55 * t64;
t46 = 0.1e1 / t47 ^ 2;
t68 = cos(qJ(1));
t80 = t46 * t68 ^ 2;
t65 = pkin(10) + qJ(5);
t62 = cos(t65);
t70 = t68 * t62;
t61 = sin(t65);
t74 = t67 * t61;
t53 = t64 * t70 + t74;
t51 = 0.1e1 / t53 ^ 2;
t71 = t68 * t61;
t73 = t67 * t62;
t52 = t64 * t71 - t73;
t79 = t51 * t52;
t78 = t54 * t64;
t58 = t63 ^ 2;
t76 = t58 / t64 ^ 2;
t57 = 0.1e1 / (t67 ^ 2 * t76 + 0.1e1);
t77 = t57 * t67;
t75 = t63 * t68;
t69 = t52 ^ 2 * t51 + 0.1e1;
t59 = 0.1e1 / t64;
t50 = 0.1e1 / t53;
t49 = (0.1e1 + t76) * t77;
t48 = 0.1e1 / t69;
t45 = 0.1e1 / t47;
t44 = 0.1e1 / (t58 * t80 + 0.1e1);
t43 = (-t50 * t61 + t62 * t79) * t48 * t75;
t42 = (t64 * t45 - (-t67 * t78 + t55 * t63 + (-t55 * t72 + t78) * t49) * t63 * t46) * t68 * t44;
t1 = [t59 * t57 * t75, t49, t49, 0, 0, 0; (-t45 * t72 - (-t55 * t58 * t59 * t77 + (t57 - 0.1e1) * t63 * t54) * t63 * t80) * t44, t42, t42, 0, 0, 0; ((-t64 * t74 - t70) * t50 - (-t64 * t73 + t71) * t79) * t48, t43, t43, 0, t69 * t48, 0;];
Ja_rot  = t1;
