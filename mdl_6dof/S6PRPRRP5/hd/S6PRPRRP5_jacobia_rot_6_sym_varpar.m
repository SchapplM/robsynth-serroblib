% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRP5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:52:31
% EndTime: 2019-02-26 19:52:31
% DurationCPUTime: 0.11s
% Computational Cost: add. (334->29), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->51)
t70 = sin(pkin(10));
t72 = cos(pkin(10));
t76 = sin(qJ(2));
t73 = cos(pkin(6));
t79 = cos(qJ(2));
t82 = t73 * t79;
t64 = t70 * t76 - t72 * t82;
t78 = cos(qJ(4));
t71 = sin(pkin(6));
t75 = sin(qJ(4));
t87 = t71 * t75;
t59 = t64 * t78 + t72 * t87;
t84 = t71 * t79;
t68 = t73 * t75 + t78 * t84;
t56 = atan2(t59, t68);
t53 = sin(t56);
t54 = cos(t56);
t47 = t53 * t59 + t54 * t68;
t46 = 0.1e1 / t47 ^ 2;
t66 = t70 * t82 + t72 * t76;
t57 = -t66 * t78 + t70 * t87;
t92 = t46 * t57;
t85 = t71 * t78;
t58 = t66 * t75 + t70 * t85;
t77 = cos(qJ(5));
t83 = t73 * t76;
t67 = -t70 * t83 + t72 * t79;
t74 = sin(qJ(5));
t89 = t67 * t74;
t52 = t58 * t77 + t89;
t50 = 0.1e1 / t52 ^ 2;
t88 = t67 * t77;
t51 = t58 * t74 - t88;
t91 = t50 * t51;
t63 = 0.1e1 / t68 ^ 2;
t90 = t59 * t63;
t86 = t71 * t76;
t81 = t51 ^ 2 * t50 + 0.1e1;
t80 = -t53 * t68 + t54 * t59;
t69 = t73 * t78 - t75 * t84;
t65 = t70 * t79 + t72 * t83;
t62 = 0.1e1 / t68;
t60 = -t64 * t75 + t72 * t85;
t55 = 0.1e1 / (t59 ^ 2 * t63 + 0.1e1);
t49 = 0.1e1 / t52;
t48 = 0.1e1 / t81;
t45 = 0.1e1 / t47;
t44 = 0.1e1 / (t57 ^ 2 * t46 + 0.1e1);
t43 = (t62 * t65 + t86 * t90) * t78 * t55;
t42 = (t60 * t62 - t69 * t90) * t55;
t1 = [0, t43, 0, t42, 0, 0; 0 (-t67 * t78 * t45 - ((t53 * t65 - t54 * t86) * t78 + t80 * t43) * t92) * t44, 0 (t58 * t45 - (t80 * t42 + t53 * t60 + t54 * t69) * t92) * t44, 0, 0; 0 ((t66 * t77 + t75 * t89) * t49 - (-t66 * t74 + t75 * t88) * t91) * t48, 0 (-t49 * t74 + t77 * t91) * t57 * t48, t81 * t48, 0;];
Ja_rot  = t1;
