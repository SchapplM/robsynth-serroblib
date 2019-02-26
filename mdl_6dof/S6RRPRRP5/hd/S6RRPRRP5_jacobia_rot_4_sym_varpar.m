% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP5_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:20
% EndTime: 2019-02-26 21:48:21
% DurationCPUTime: 0.15s
% Computational Cost: add. (374->26), mult. (1051->69), div. (55->9), fcn. (1490->13), ass. (0->45)
t73 = sin(qJ(1));
t76 = cos(qJ(1));
t68 = sin(pkin(11));
t70 = cos(pkin(11));
t75 = cos(qJ(2));
t83 = cos(pkin(6));
t80 = t75 * t83;
t72 = sin(qJ(2));
t81 = t72 * t83;
t77 = -t68 * t81 + t70 * t80;
t78 = t75 * t68 + t72 * t70;
t54 = -t73 * t78 + t76 * t77;
t65 = t72 * t68 - t75 * t70;
t69 = sin(pkin(6));
t62 = t65 * t69;
t48 = atan2(t54, t62);
t46 = cos(t48);
t88 = t46 * t54;
t64 = t68 * t80 + t70 * t81;
t58 = -t73 * t64 - t76 * t65;
t71 = sin(qJ(4));
t74 = cos(qJ(4));
t85 = t69 * t73;
t52 = t58 * t74 + t71 * t85;
t50 = 0.1e1 / t52 ^ 2;
t51 = t58 * t71 - t74 * t85;
t87 = t50 * t51;
t45 = sin(t48);
t43 = t45 * t54 + t46 * t62;
t42 = 0.1e1 / t43 ^ 2;
t56 = -t73 * t77 - t76 * t78;
t86 = t56 ^ 2 * t42;
t84 = t69 * t76;
t82 = t51 ^ 2 * t50 + 0.1e1;
t79 = -t76 * t64 + t73 * t65;
t63 = t78 * t69;
t61 = 0.1e1 / t62 ^ 2;
t60 = 0.1e1 / t62;
t49 = 0.1e1 / t52;
t47 = 0.1e1 / (t54 ^ 2 * t61 + 0.1e1);
t44 = 0.1e1 / t82;
t41 = 0.1e1 / t43;
t40 = 0.1e1 / (0.1e1 + t86);
t39 = (-t54 * t61 * t63 + t60 * t79) * t47;
t1 = [t56 * t60 * t47, t39, 0, 0, 0, 0; (t54 * t41 + (t45 + (t60 * t88 - t45) * t47) * t86) * t40 (t58 * t41 + (t45 * t79 + t46 * t63 + (-t45 * t62 + t88) * t39) * t56 * t42) * t40, 0, 0, 0, 0; ((t71 * t79 - t74 * t84) * t49 - (t71 * t84 + t74 * t79) * t87) * t44 (t71 * t49 - t74 * t87) * t56 * t44, 0, t82 * t44, 0, 0;];
Ja_rot  = t1;
