% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR15_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobia_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobia_rot_3_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:47
% EndTime: 2019-02-26 22:38:47
% DurationCPUTime: 0.14s
% Computational Cost: add. (304->29), mult. (885->75), div. (55->9), fcn. (1254->13), ass. (0->49)
t69 = cos(pkin(6));
t74 = cos(qJ(2));
t75 = cos(qJ(1));
t79 = t74 * t75;
t71 = sin(qJ(2));
t72 = sin(qJ(1));
t81 = t72 * t71;
t60 = -t69 * t79 + t81;
t66 = sin(pkin(7));
t68 = cos(pkin(7));
t67 = sin(pkin(6));
t83 = t67 * t75;
t54 = -t60 * t66 + t68 * t83;
t59 = -t66 * t67 * t74 + t68 * t69;
t53 = atan2(t54, t59);
t50 = sin(t53);
t51 = cos(t53);
t44 = t50 * t54 + t51 * t59;
t43 = 0.1e1 / t44 ^ 2;
t80 = t72 * t74;
t82 = t71 * t75;
t62 = -t69 * t80 - t82;
t84 = t67 * t72;
t55 = t62 * t66 - t68 * t84;
t90 = t43 * t55 ^ 2;
t70 = sin(qJ(3));
t76 = t62 * t68 + t66 * t84;
t63 = -t69 * t81 + t79;
t73 = cos(qJ(3));
t86 = t63 * t73;
t49 = t76 * t70 + t86;
t47 = 0.1e1 / t49 ^ 2;
t87 = t63 * t70;
t48 = -t76 * t73 + t87;
t89 = t47 * t48;
t88 = t51 * t54;
t85 = t67 * t71;
t78 = t47 * t48 ^ 2 + 0.1e1;
t77 = t60 * t68 + t66 * t83;
t61 = -t69 * t82 - t80;
t58 = 0.1e1 / t59 ^ 2;
t57 = 0.1e1 / t59;
t52 = 0.1e1 / (t54 ^ 2 * t58 + 0.1e1);
t46 = 0.1e1 / t49;
t45 = 0.1e1 / t78;
t42 = 0.1e1 / t44;
t41 = 0.1e1 / (0.1e1 + t90);
t40 = (-t54 * t58 * t85 + t57 * t61) * t66 * t52;
t1 = [t55 * t57 * t52, t40, 0, 0, 0, 0; (t54 * t42 + (t50 + (t57 * t88 - t50) * t52) * t90) * t41 (t63 * t66 * t42 + ((t50 * t61 + t51 * t85) * t66 + (-t50 * t59 + t88) * t40) * t55 * t43) * t41, 0, 0, 0, 0; ((t61 * t70 - t77 * t73) * t46 - (t61 * t73 + t77 * t70) * t89) * t45 ((t62 * t70 + t68 * t86) * t46 - (t62 * t73 - t68 * t87) * t89) * t45, t78 * t45, 0, 0, 0;];
Ja_rot  = t1;
