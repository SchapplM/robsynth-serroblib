% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:32:34
% EndTime: 2019-02-26 21:32:34
% DurationCPUTime: 0.09s
% Computational Cost: add. (199->25), mult. (409->67), div. (57->9), fcn. (587->11), ass. (0->41)
t69 = cos(pkin(10));
t71 = sin(qJ(1));
t72 = cos(qJ(2));
t68 = sin(pkin(10));
t73 = cos(qJ(1));
t76 = t73 * t68;
t56 = -t71 * t69 + t72 * t76;
t75 = t73 * t69;
t57 = t71 * t68 + t72 * t75;
t67 = qJ(5) + qJ(6);
t62 = sin(t67);
t63 = cos(t67);
t48 = t56 * t62 + t57 * t63;
t46 = 0.1e1 / t48 ^ 2;
t47 = -t56 * t63 + t57 * t62;
t83 = t46 * t47;
t70 = sin(qJ(2));
t78 = t71 * t70;
t61 = atan2(t78, t72);
t58 = sin(t61);
t59 = cos(t61);
t52 = t58 * t78 + t59 * t72;
t51 = 0.1e1 / t52 ^ 2;
t82 = t51 * t73 ^ 2;
t64 = t70 ^ 2;
t81 = t64 / t72 ^ 2;
t80 = t70 * t73;
t60 = 0.1e1 / (t71 ^ 2 * t81 + 0.1e1);
t79 = t71 * t60;
t77 = t71 * t72;
t74 = t47 ^ 2 * t46 + 0.1e1;
t65 = 0.1e1 / t72;
t55 = -t69 * t77 + t76;
t54 = -t68 * t77 - t75;
t53 = (0.1e1 + t81) * t79;
t50 = 0.1e1 / t52;
t49 = 0.1e1 / (t64 * t82 + 0.1e1);
t45 = 0.1e1 / t48;
t44 = 0.1e1 / t74;
t43 = t74 * t44;
t1 = [t65 * t60 * t80, t53, 0, 0, 0, 0; (t50 * t78 + (t59 * t64 * t65 * t79 + (-t60 + 0.1e1) * t70 * t58) * t70 * t82) * t49 (-t72 * t50 + (t58 * t77 - t59 * t70 + (-t58 * t72 + t59 * t78) * t53) * t70 * t51) * t73 * t49, 0, 0, 0, 0; ((-t54 * t63 + t55 * t62) * t45 - (t54 * t62 + t55 * t63) * t83) * t44 ((-t62 * t69 + t63 * t68) * t45 - (-t62 * t68 - t63 * t69) * t83) * t44 * t80, 0, 0, t43, t43;];
Ja_rot  = t1;
