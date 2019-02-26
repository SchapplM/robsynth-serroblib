% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR11_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:48
% EndTime: 2019-02-26 21:59:48
% DurationCPUTime: 0.11s
% Computational Cost: add. (259->21), mult. (278->54), div. (62->9), fcn. (402->9), ass. (0->37)
t60 = sin(qJ(2));
t61 = sin(qJ(1));
t62 = cos(qJ(2));
t68 = t61 * t62;
t52 = atan2(-t68, t60);
t50 = sin(t52);
t51 = cos(t52);
t44 = -t50 * t68 + t51 * t60;
t43 = 0.1e1 / t44 ^ 2;
t63 = cos(qJ(1));
t75 = t43 * t63 ^ 2;
t56 = qJ(4) + qJ(5) + qJ(6);
t54 = sin(t56);
t66 = t63 * t54;
t55 = cos(t56);
t69 = t61 * t55;
t49 = t60 * t66 + t69;
t47 = 0.1e1 / t49 ^ 2;
t65 = t63 * t55;
t70 = t61 * t54;
t48 = -t60 * t65 + t70;
t74 = t47 * t48;
t73 = t50 * t60;
t59 = t62 ^ 2;
t72 = 0.1e1 / t60 ^ 2 * t59;
t53 = 0.1e1 / (t61 ^ 2 * t72 + 0.1e1);
t71 = t61 * t53;
t67 = t62 * t63;
t64 = t48 ^ 2 * t47 + 0.1e1;
t57 = 0.1e1 / t60;
t46 = 0.1e1 / t49;
t45 = (0.1e1 + t72) * t71;
t42 = 0.1e1 / t44;
t41 = 0.1e1 / (t59 * t75 + 0.1e1);
t40 = 0.1e1 / t64;
t39 = t64 * t40;
t1 = [-t57 * t53 * t67, t45, 0, 0, 0, 0; (-t42 * t68 - (t51 * t57 * t59 * t71 + (t53 - 0.1e1) * t62 * t50) * t62 * t75) * t41 (-t60 * t42 - (t61 * t73 + t51 * t62 + (-t51 * t68 - t73) * t45) * t62 * t43) * t63 * t41, 0, 0, 0, 0; ((t60 * t69 + t66) * t46 - (-t60 * t70 + t65) * t74) * t40 (-t46 * t55 - t54 * t74) * t40 * t67, 0, t39, t39, t39;];
Ja_rot  = t1;
