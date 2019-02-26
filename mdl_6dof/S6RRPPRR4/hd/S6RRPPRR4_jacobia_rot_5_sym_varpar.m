% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:11
% EndTime: 2019-02-26 21:30:11
% DurationCPUTime: 0.18s
% Computational Cost: add. (374->26), mult. (1051->65), div. (55->9), fcn. (1490->13), ass. (0->43)
t72 = sin(pkin(11));
t74 = cos(pkin(11));
t77 = sin(qJ(2));
t80 = cos(qJ(2));
t69 = t77 * t72 - t80 * t74;
t78 = sin(qJ(1));
t81 = cos(qJ(1));
t75 = cos(pkin(6));
t84 = t80 * t72 + t77 * t74;
t83 = t84 * t75;
t62 = -t81 * t69 - t78 * t83;
t58 = -t78 * t69 + t81 * t83;
t73 = sin(pkin(6));
t68 = t84 * t73;
t52 = atan2(-t58, t68);
t50 = cos(t52);
t91 = t50 * t58;
t82 = t69 * t75;
t61 = t78 * t82 - t81 * t84;
t76 = sin(qJ(5));
t79 = cos(qJ(5));
t88 = t73 * t78;
t56 = -t61 * t76 + t79 * t88;
t54 = 0.1e1 / t56 ^ 2;
t55 = t61 * t79 + t76 * t88;
t90 = t54 * t55;
t49 = sin(t52);
t47 = -t49 * t58 + t50 * t68;
t46 = 0.1e1 / t47 ^ 2;
t89 = t62 ^ 2 * t46;
t87 = t73 * t81;
t85 = t55 ^ 2 * t54 + 0.1e1;
t67 = t69 * t73;
t66 = 0.1e1 / t68 ^ 2;
t65 = 0.1e1 / t68;
t57 = -t78 * t84 - t81 * t82;
t53 = 0.1e1 / t56;
t51 = 0.1e1 / (t58 ^ 2 * t66 + 0.1e1);
t48 = 0.1e1 / t85;
t45 = 0.1e1 / t47;
t44 = 0.1e1 / (0.1e1 + t89);
t43 = (-t58 * t66 * t67 - t57 * t65) * t51;
t1 = [-t62 * t65 * t51, t43, 0, 0, 0, 0; (-t58 * t45 - (-t49 + (t65 * t91 + t49) * t51) * t89) * t44 (t61 * t45 - (-t49 * t57 - t50 * t67 + (-t49 * t68 - t91) * t43) * t62 * t46) * t44, 0, 0, 0, 0; ((-t57 * t79 + t76 * t87) * t53 - (t57 * t76 + t79 * t87) * t90) * t48 (-t79 * t53 - t76 * t90) * t62 * t48, 0, 0, t85 * t48, 0;];
Ja_rot  = t1;
