% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRRRR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobia_rot_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobia_rot_5_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:29:14
% EndTime: 2019-07-18 13:29:14
% DurationCPUTime: 0.08s
% Computational Cost: add. (323->20), mult. (306->52), div. (109->11), fcn. (491->9), ass. (0->40)
t59 = qJ(3) + qJ(4);
t56 = cos(t59);
t55 = sin(t59);
t61 = sin(qJ(2));
t71 = t61 * t55;
t50 = atan2(-t56, t71);
t49 = cos(t50);
t48 = sin(t50);
t75 = t48 * t56;
t41 = t49 * t71 - t75;
t40 = 0.1e1 / t41 ^ 2;
t63 = cos(qJ(2));
t77 = t40 * t63 ^ 2;
t62 = cos(qJ(5));
t67 = t63 * t62;
t60 = sin(qJ(5));
t70 = t61 * t60;
t47 = t56 * t67 + t70;
t45 = 0.1e1 / t47 ^ 2;
t68 = t63 * t60;
t69 = t61 * t62;
t46 = t56 * t68 - t69;
t76 = t45 * t46;
t74 = t49 * t56;
t54 = t56 ^ 2;
t64 = t55 ^ 2;
t73 = 0.1e1 / t64 * t54;
t58 = 0.1e1 / t61 ^ 2;
t51 = 0.1e1 / (t58 * t73 + 0.1e1);
t72 = 0.1e1 / t61 * t51;
t66 = t51 / t55 * t58;
t65 = t46 ^ 2 * t45 + 0.1e1;
t44 = 0.1e1 / t47;
t43 = 0.1e1 / t65;
t42 = (0.1e1 + t73) * t72;
t39 = 0.1e1 / t41;
t38 = 0.1e1 / (t64 * t77 + 0.1e1);
t37 = (-t44 * t60 + t62 * t76) * t63 * t55 * t43;
t36 = (t56 * t39 - (t61 * t74 + t48 * t55 + (-t48 * t71 - t74) * t42) * t55 * t40) * t63 * t38;
t1 = [0, t63 * t56 * t66, t42, t42, 0; 0, (-t39 * t71 - (-t72 * t75 + (-t54 * t66 + t55) * t49) * t55 * t77) * t38, t36, t36, 0; 0, ((-t56 * t70 - t67) * t44 - (-t56 * t69 + t68) * t76) * t43, t37, t37, t65 * t43;];
Ja_rot  = t1;
