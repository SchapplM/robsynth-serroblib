% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR10_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:48
% EndTime: 2019-02-26 22:35:48
% DurationCPUTime: 0.09s
% Computational Cost: add. (252->25), mult. (499->64), div. (77->11), fcn. (736->11), ass. (0->42)
t69 = sin(qJ(2));
t71 = cos(qJ(2));
t72 = cos(qJ(1));
t70 = sin(qJ(1));
t76 = cos(pkin(6));
t74 = t70 * t76;
t60 = -t69 * t74 + t72 * t71;
t67 = qJ(3) + qJ(4);
t62 = sin(t67);
t63 = cos(t67);
t68 = sin(pkin(6));
t79 = t68 * t70;
t51 = t60 * t63 + t62 * t79;
t49 = 0.1e1 / t51 ^ 2;
t50 = t60 * t62 - t63 * t79;
t83 = t49 * t50;
t73 = t72 * t76;
t56 = t70 * t69 - t71 * t73;
t78 = t68 * t71;
t54 = atan2(-t56, -t78);
t53 = cos(t54);
t82 = t53 * t56;
t52 = sin(t54);
t46 = -t52 * t56 - t53 * t78;
t45 = 0.1e1 / t46 ^ 2;
t59 = t72 * t69 + t71 * t74;
t81 = t59 ^ 2 * t45;
t64 = 0.1e1 / t68;
t65 = 0.1e1 / t71;
t80 = t64 * t65;
t77 = t68 * t72;
t75 = t50 ^ 2 * t49 + 0.1e1;
t66 = 0.1e1 / t71 ^ 2;
t58 = t69 * t73 + t70 * t71;
t55 = 0.1e1 / (0.1e1 + t56 ^ 2 / t68 ^ 2 * t66);
t48 = 0.1e1 / t51;
t47 = 0.1e1 / t75;
t44 = 0.1e1 / t46;
t43 = 0.1e1 / (0.1e1 + t81);
t42 = (t56 * t66 * t69 + t58 * t65) * t64 * t55;
t41 = t75 * t47;
t1 = [t59 * t55 * t80, t42, 0, 0, 0, 0; (-t56 * t44 - (-t52 + (-t80 * t82 + t52) * t55) * t81) * t43 (t60 * t44 - (t53 * t68 * t69 - t52 * t58 + (t52 * t78 - t82) * t42) * t59 * t45) * t43, 0, 0, 0, 0; ((-t58 * t62 - t63 * t77) * t48 - (-t58 * t63 + t62 * t77) * t83) * t47 (-t62 * t48 + t63 * t83) * t59 * t47, t41, t41, 0, 0;];
Ja_rot  = t1;
