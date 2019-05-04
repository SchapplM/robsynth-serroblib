% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR14V3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobia_rot_5_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:04
% EndTime: 2019-04-12 15:12:04
% DurationCPUTime: 0.15s
% Computational Cost: add. (216->32), mult. (662->85), div. (108->11), fcn. (985->11), ass. (0->45)
t59 = sin(qJ(4));
t63 = cos(qJ(4));
t65 = cos(qJ(1));
t67 = t65 * t63;
t61 = sin(qJ(1));
t64 = cos(qJ(2));
t69 = t61 * t64;
t47 = t59 * t69 + t67;
t60 = sin(qJ(2));
t73 = t60 * t59;
t46 = atan2(-t47, t73);
t43 = sin(t46);
t44 = cos(t46);
t37 = -t43 * t47 + t44 * t73;
t36 = 0.1e1 / t37 ^ 2;
t68 = t65 * t59;
t50 = -t61 * t63 + t64 * t68;
t78 = t36 * t50;
t51 = t61 * t59 + t64 * t67;
t58 = sin(qJ(5));
t62 = cos(qJ(5));
t70 = t60 * t65;
t42 = t51 * t62 + t58 * t70;
t40 = 0.1e1 / t42 ^ 2;
t41 = t51 * t58 - t62 * t70;
t77 = t40 * t41;
t76 = t44 * t47;
t75 = t50 ^ 2 * t36;
t54 = 0.1e1 / t59;
t56 = 0.1e1 / t60;
t74 = t54 * t56;
t72 = t60 * t61;
t71 = t60 * t63;
t66 = t41 ^ 2 * t40 + 0.1e1;
t57 = 0.1e1 / t60 ^ 2;
t55 = 0.1e1 / t59 ^ 2;
t49 = t63 * t69 - t68;
t45 = 0.1e1 / (t47 ^ 2 * t57 * t55 + 0.1e1);
t39 = 0.1e1 / t42;
t38 = 0.1e1 / t66;
t35 = 0.1e1 / t37;
t34 = (t47 * t54 * t57 * t64 + t61) * t45;
t33 = 0.1e1 / (0.1e1 + t75);
t32 = (t47 * t55 * t63 - t49 * t54) * t56 * t45;
t1 = [-t50 * t45 * t74, t34, 0, t32, 0, 0; (-t47 * t35 - (-t43 + (t74 * t76 + t43) * t45) * t75) * t33 (t34 * t76 * t78 + (-t35 * t70 - (t44 * t64 + (-t34 * t60 + t72) * t43) * t78) * t59) * t33, 0 (t51 * t35 - (t44 * t71 - t43 * t49 + (-t43 * t73 - t76) * t32) * t78) * t33, 0, 0; ((-t49 * t58 + t62 * t72) * t39 - (-t49 * t62 - t58 * t72) * t77) * t38 ((-t58 * t71 - t62 * t64) * t39 - (t58 * t64 - t62 * t71) * t77) * t38 * t65, 0 (-t58 * t39 + t62 * t77) * t50 * t38, t66 * t38, 0;];
Ja_rot  = t1;
