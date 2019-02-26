% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR9_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:58:50
% EndTime: 2019-02-26 21:58:50
% DurationCPUTime: 0.08s
% Computational Cost: add. (314->25), mult. (499->64), div. (77->11), fcn. (736->11), ass. (0->42)
t71 = sin(qJ(2));
t73 = cos(qJ(2));
t74 = cos(qJ(1));
t72 = sin(qJ(1));
t78 = cos(pkin(6));
t76 = t72 * t78;
t62 = -t71 * t76 + t74 * t73;
t66 = pkin(12) + qJ(4) + qJ(5);
t64 = sin(t66);
t65 = cos(t66);
t70 = sin(pkin(6));
t81 = t70 * t72;
t53 = t62 * t65 + t64 * t81;
t51 = 0.1e1 / t53 ^ 2;
t52 = t62 * t64 - t65 * t81;
t85 = t51 * t52;
t75 = t74 * t78;
t58 = t72 * t71 - t73 * t75;
t80 = t70 * t73;
t56 = atan2(-t58, -t80);
t55 = cos(t56);
t84 = t55 * t58;
t54 = sin(t56);
t48 = -t54 * t58 - t55 * t80;
t47 = 0.1e1 / t48 ^ 2;
t61 = t74 * t71 + t73 * t76;
t83 = t61 ^ 2 * t47;
t67 = 0.1e1 / t70;
t68 = 0.1e1 / t73;
t82 = t67 * t68;
t79 = t70 * t74;
t77 = t52 ^ 2 * t51 + 0.1e1;
t69 = 0.1e1 / t73 ^ 2;
t60 = t71 * t75 + t72 * t73;
t57 = 0.1e1 / (0.1e1 + t58 ^ 2 / t70 ^ 2 * t69);
t50 = 0.1e1 / t53;
t49 = 0.1e1 / t77;
t46 = 0.1e1 / t48;
t45 = 0.1e1 / (0.1e1 + t83);
t44 = (t58 * t69 * t71 + t60 * t68) * t67 * t57;
t43 = t77 * t49;
t1 = [t61 * t57 * t82, t44, 0, 0, 0, 0; (-t58 * t46 - (-t54 + (-t82 * t84 + t54) * t57) * t83) * t45 (t62 * t46 - (t55 * t70 * t71 - t54 * t60 + (t54 * t80 - t84) * t44) * t61 * t47) * t45, 0, 0, 0, 0; ((-t60 * t64 - t65 * t79) * t50 - (-t60 * t65 + t64 * t79) * t85) * t49 (-t64 * t50 + t65 * t85) * t61 * t49, 0, t43, t43, 0;];
Ja_rot  = t1;
