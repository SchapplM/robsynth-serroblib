% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:12
% EndTime: 2019-02-26 19:54:12
% DurationCPUTime: 0.22s
% Computational Cost: add. (632->32), mult. (1727->84), div. (65->9), fcn. (2425->15), ass. (0->54)
t85 = sin(pkin(6));
t93 = cos(qJ(4));
t100 = t85 * t93;
t88 = cos(pkin(6));
t83 = sin(pkin(12));
t86 = cos(pkin(12));
t91 = sin(qJ(2));
t94 = cos(qJ(2));
t96 = t94 * t83 + t91 * t86;
t78 = t96 * t88;
t79 = t91 * t83 - t94 * t86;
t84 = sin(pkin(11));
t87 = cos(pkin(11));
t67 = t87 * t78 - t84 * t79;
t90 = sin(qJ(4));
t61 = t87 * t100 + t67 * t90;
t77 = t96 * t85;
t73 = t77 * t90 - t88 * t93;
t60 = atan2(-t61, t73);
t57 = sin(t60);
t58 = cos(t60);
t51 = -t57 * t61 + t58 * t73;
t50 = 0.1e1 / t51 ^ 2;
t97 = -t84 * t78 - t87 * t79;
t64 = -t84 * t100 + t90 * t97;
t105 = t50 * t64;
t101 = t85 * t90;
t65 = t84 * t101 + t93 * t97;
t95 = t79 * t88;
t69 = t84 * t95 - t87 * t96;
t89 = sin(qJ(5));
t92 = cos(qJ(5));
t56 = t65 * t92 - t69 * t89;
t54 = 0.1e1 / t56 ^ 2;
t55 = t65 * t89 + t69 * t92;
t104 = t54 * t55;
t72 = 0.1e1 / t73 ^ 2;
t103 = t61 * t72;
t102 = t69 * t93;
t99 = t55 ^ 2 * t54 + 0.1e1;
t98 = -t57 * t73 - t58 * t61;
t76 = t79 * t85;
t74 = t77 * t93 + t88 * t90;
t71 = 0.1e1 / t73;
t66 = -t84 * t96 - t87 * t95;
t63 = -t87 * t101 + t67 * t93;
t59 = 0.1e1 / (t61 ^ 2 * t72 + 0.1e1);
t53 = 0.1e1 / t56;
t52 = 0.1e1 / t99;
t49 = 0.1e1 / t51;
t48 = 0.1e1 / (t64 ^ 2 * t50 + 0.1e1);
t47 = (-t76 * t103 - t66 * t71) * t90 * t59;
t46 = (t74 * t103 - t63 * t71) * t59;
t1 = [0, t47, 0, t46, 0, 0; 0 (t69 * t90 * t49 - ((-t57 * t66 - t58 * t76) * t90 + t98 * t47) * t105) * t48, 0 (t65 * t49 - (t98 * t46 - t57 * t63 + t58 * t74) * t105) * t48, 0, 0; 0 ((t89 * t102 - t92 * t97) * t53 - (t92 * t102 + t89 * t97) * t104) * t52, 0 (t92 * t104 - t53 * t89) * t64 * t52, t99 * t52, 0;];
Ja_rot  = t1;
