% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:44
% EndTime: 2019-02-26 19:47:44
% DurationCPUTime: 0.15s
% Computational Cost: add. (481->27), mult. (1334->73), div. (57->9), fcn. (1884->13), ass. (0->47)
t78 = cos(pkin(6));
t73 = sin(pkin(11));
t76 = cos(pkin(11));
t80 = sin(qJ(2));
t82 = cos(qJ(2));
t84 = t82 * t73 + t80 * t76;
t69 = t84 * t78;
t70 = t80 * t73 - t82 * t76;
t74 = sin(pkin(10));
t77 = cos(pkin(10));
t58 = t77 * t69 - t74 * t70;
t79 = sin(qJ(4));
t75 = sin(pkin(6));
t81 = cos(qJ(4));
t86 = t75 * t81;
t50 = t58 * t79 + t77 * t86;
t68 = t84 * t75;
t64 = t68 * t79 - t78 * t81;
t49 = atan2(-t50, t64);
t46 = sin(t49);
t47 = cos(t49);
t44 = -t46 * t50 + t47 * t64;
t43 = 0.1e1 / t44 ^ 2;
t61 = -t74 * t69 - t77 * t70;
t53 = t61 * t79 - t74 * t86;
t89 = t43 * t53;
t63 = 0.1e1 / t64 ^ 2;
t88 = t50 * t63;
t87 = t75 * t79;
t85 = -t46 * t64 - t47 * t50;
t83 = t70 * t78;
t67 = t70 * t75;
t65 = t68 * t81 + t78 * t79;
t62 = 0.1e1 / t64;
t59 = t74 * t83 - t77 * t84;
t57 = -t74 * t84 - t77 * t83;
t56 = 0.1e1 / t59 ^ 2;
t55 = 0.1e1 / t59;
t54 = t61 * t81 + t74 * t87;
t52 = t58 * t81 - t77 * t87;
t48 = 0.1e1 / (t50 ^ 2 * t63 + 0.1e1);
t45 = 0.1e1 / (t54 ^ 2 * t56 + 0.1e1);
t42 = 0.1e1 / t44;
t41 = 0.1e1 / (t53 ^ 2 * t43 + 0.1e1);
t40 = (-t57 * t62 - t67 * t88) * t79 * t48;
t39 = (-t52 * t62 + t65 * t88) * t48;
t1 = [0, t40, 0, t39, 0, 0; 0 (t59 * t79 * t42 - ((-t46 * t57 - t47 * t67) * t79 + t85 * t40) * t89) * t41, 0 (t54 * t42 - (t85 * t39 - t46 * t52 + t47 * t65) * t89) * t41, 0, 0; 0 (-t61 * t54 * t56 - t59 * t81 * t55) * t45, 0, t53 * t55 * t45, 0, 0;];
Ja_rot  = t1;
