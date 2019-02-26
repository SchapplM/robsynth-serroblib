% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR12_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:00:24
% EndTime: 2019-02-26 22:00:24
% DurationCPUTime: 0.09s
% Computational Cost: add. (228->25), mult. (499->64), div. (77->11), fcn. (736->11), ass. (0->42)
t72 = sin(qJ(2));
t74 = cos(qJ(2));
t75 = cos(qJ(1));
t73 = sin(qJ(1));
t79 = cos(pkin(6));
t77 = t73 * t79;
t61 = t75 * t72 + t74 * t77;
t70 = qJ(4) + qJ(5);
t65 = sin(t70);
t66 = cos(t70);
t71 = sin(pkin(6));
t81 = t71 * t73;
t53 = t61 * t65 + t66 * t81;
t51 = 0.1e1 / t53 ^ 2;
t52 = -t61 * t66 + t65 * t81;
t86 = t51 * t52;
t76 = t75 * t79;
t59 = t72 * t76 + t73 * t74;
t82 = t71 * t72;
t57 = atan2(-t59, t82);
t55 = cos(t57);
t85 = t55 * t59;
t54 = sin(t57);
t48 = -t54 * t59 + t55 * t82;
t47 = 0.1e1 / t48 ^ 2;
t62 = -t72 * t77 + t75 * t74;
t84 = t62 ^ 2 * t47;
t67 = 0.1e1 / t71;
t68 = 0.1e1 / t72;
t83 = t67 * t68;
t80 = t71 * t75;
t78 = t52 ^ 2 * t51 + 0.1e1;
t69 = 0.1e1 / t72 ^ 2;
t58 = t73 * t72 - t74 * t76;
t56 = 0.1e1 / (0.1e1 + t59 ^ 2 / t71 ^ 2 * t69);
t50 = 0.1e1 / t53;
t49 = 0.1e1 / t78;
t46 = 0.1e1 / t48;
t45 = 0.1e1 / (0.1e1 + t84);
t44 = (t59 * t69 * t74 + t58 * t68) * t67 * t56;
t43 = t78 * t49;
t1 = [-t62 * t56 * t83, t44, 0, 0, 0, 0; (-t59 * t46 - (-t54 + (t83 * t85 + t54) * t56) * t84) * t45 (-t61 * t46 - (t55 * t71 * t74 + t54 * t58 + (-t54 * t82 - t85) * t44) * t62 * t47) * t45, 0, 0, 0, 0; ((t58 * t66 + t65 * t80) * t50 - (-t58 * t65 + t66 * t80) * t86) * t49 (-t66 * t50 - t65 * t86) * t62 * t49, 0, t43, t43, 0;];
Ja_rot  = t1;
