% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:16
% EndTime: 2019-02-26 19:48:16
% DurationCPUTime: 0.13s
% Computational Cost: add. (596->31), mult. (866->78), div. (60->9), fcn. (1237->13), ass. (0->51)
t70 = sin(pkin(10));
t73 = cos(pkin(10));
t76 = cos(qJ(2));
t74 = cos(pkin(6));
t75 = sin(qJ(2));
t79 = t74 * t75;
t62 = t70 * t76 + t73 * t79;
t68 = pkin(11) + qJ(4);
t66 = sin(t68);
t67 = cos(t68);
t71 = sin(pkin(6));
t82 = t71 * t73;
t52 = t62 * t66 + t67 * t82;
t81 = t71 * t75;
t59 = t66 * t81 - t74 * t67;
t51 = atan2(-t52, t59);
t48 = sin(t51);
t49 = cos(t51);
t42 = -t48 * t52 + t49 * t59;
t41 = 0.1e1 / t42 ^ 2;
t64 = -t70 * t79 + t73 * t76;
t83 = t70 * t71;
t55 = t64 * t66 - t67 * t83;
t88 = t41 * t55;
t56 = t64 * t67 + t66 * t83;
t72 = cos(pkin(12));
t78 = t74 * t76;
t63 = t70 * t78 + t73 * t75;
t69 = sin(pkin(12));
t85 = t63 * t69;
t47 = t56 * t72 + t85;
t45 = 0.1e1 / t47 ^ 2;
t84 = t63 * t72;
t46 = t56 * t69 - t84;
t87 = t45 * t46;
t58 = 0.1e1 / t59 ^ 2;
t86 = t52 * t58;
t80 = t71 * t76;
t77 = -t48 * t59 - t49 * t52;
t61 = -t70 * t75 + t73 * t78;
t60 = t74 * t66 + t67 * t81;
t57 = 0.1e1 / t59;
t54 = t62 * t67 - t66 * t82;
t50 = 0.1e1 / (t52 ^ 2 * t58 + 0.1e1);
t44 = 0.1e1 / t47;
t43 = 0.1e1 / (t46 ^ 2 * t45 + 0.1e1);
t40 = 0.1e1 / t42;
t39 = 0.1e1 / (t55 ^ 2 * t41 + 0.1e1);
t38 = (-t57 * t61 + t80 * t86) * t66 * t50;
t37 = (-t54 * t57 + t60 * t86) * t50;
t1 = [0, t38, 0, t37, 0, 0; 0 (-t63 * t66 * t40 - ((-t48 * t61 + t49 * t80) * t66 + t77 * t38) * t88) * t39, 0 (t56 * t40 - (t77 * t37 - t48 * t54 + t49 * t60) * t88) * t39, 0, 0; 0 ((-t64 * t72 - t67 * t85) * t44 - (t64 * t69 - t67 * t84) * t87) * t43, 0 (-t44 * t69 + t72 * t87) * t55 * t43, 0, 0;];
Ja_rot  = t1;
