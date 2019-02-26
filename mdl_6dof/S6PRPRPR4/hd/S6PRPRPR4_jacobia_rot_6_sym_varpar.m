% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_rot = S6PRPRPR4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:15
% EndTime: 2019-02-26 19:48:15
% DurationCPUTime: 0.13s
% Computational Cost: add. (689->32), mult. (949->80), div. (65->9), fcn. (1349->13), ass. (0->52)
t77 = sin(pkin(10));
t79 = cos(pkin(10));
t82 = cos(qJ(2));
t80 = cos(pkin(6));
t81 = sin(qJ(2));
t86 = t80 * t81;
t67 = t77 * t82 + t79 * t86;
t76 = pkin(11) + qJ(4);
t72 = sin(t76);
t74 = cos(t76);
t78 = sin(pkin(6));
t89 = t78 * t79;
t57 = t67 * t72 + t74 * t89;
t88 = t78 * t81;
t64 = t72 * t88 - t80 * t74;
t56 = atan2(-t57, t64);
t53 = sin(t56);
t54 = cos(t56);
t47 = -t53 * t57 + t54 * t64;
t46 = 0.1e1 / t47 ^ 2;
t69 = -t77 * t86 + t79 * t82;
t90 = t77 * t78;
t60 = t69 * t72 - t74 * t90;
t94 = t46 * t60;
t61 = t69 * t74 + t72 * t90;
t85 = t80 * t82;
t68 = t77 * t85 + t79 * t81;
t75 = pkin(12) + qJ(6);
t71 = sin(t75);
t73 = cos(t75);
t52 = t61 * t73 + t68 * t71;
t50 = 0.1e1 / t52 ^ 2;
t51 = t61 * t71 - t68 * t73;
t93 = t50 * t51;
t63 = 0.1e1 / t64 ^ 2;
t92 = t57 * t63;
t91 = t68 * t74;
t87 = t78 * t82;
t84 = t51 ^ 2 * t50 + 0.1e1;
t83 = -t53 * t64 - t54 * t57;
t66 = -t77 * t81 + t79 * t85;
t65 = t80 * t72 + t74 * t88;
t62 = 0.1e1 / t64;
t59 = t67 * t74 - t72 * t89;
t55 = 0.1e1 / (t57 ^ 2 * t63 + 0.1e1);
t49 = 0.1e1 / t52;
t48 = 0.1e1 / t84;
t45 = 0.1e1 / t47;
t44 = 0.1e1 / (t60 ^ 2 * t46 + 0.1e1);
t43 = (-t62 * t66 + t87 * t92) * t72 * t55;
t42 = (-t59 * t62 + t65 * t92) * t55;
t1 = [0, t43, 0, t42, 0, 0; 0 (-t68 * t72 * t45 - ((-t53 * t66 + t54 * t87) * t72 + t83 * t43) * t94) * t44, 0 (t61 * t45 - (t83 * t42 - t53 * t59 + t54 * t65) * t94) * t44, 0, 0; 0 ((-t69 * t73 - t71 * t91) * t49 - (t69 * t71 - t73 * t91) * t93) * t48, 0 (-t49 * t71 + t73 * t93) * t60 * t48, 0, t84 * t48;];
Ja_rot  = t1;
