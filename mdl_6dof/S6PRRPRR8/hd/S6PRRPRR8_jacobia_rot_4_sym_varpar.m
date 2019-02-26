% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR8_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobia_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:17
% EndTime: 2019-02-26 20:08:17
% DurationCPUTime: 0.16s
% Computational Cost: add. (473->35), mult. (1448->85), div. (57->9), fcn. (1990->13), ass. (0->52)
t76 = cos(pkin(12));
t77 = cos(pkin(7));
t73 = sin(pkin(12));
t80 = sin(qJ(2));
t78 = cos(pkin(6));
t82 = cos(qJ(2));
t90 = t78 * t82;
t84 = -t73 * t80 + t76 * t90;
t74 = sin(pkin(7));
t75 = sin(pkin(6));
t93 = t75 * t74;
t98 = -t76 * t93 + t84 * t77;
t91 = t78 * t80;
t69 = t73 * t82 + t76 * t91;
t79 = sin(qJ(3));
t81 = cos(qJ(3));
t54 = t69 * t79 - t98 * t81;
t92 = t77 * t81;
t94 = t74 * t78;
t64 = -t81 * t94 + (t79 * t80 - t82 * t92) * t75;
t53 = atan2(-t54, t64);
t49 = sin(t53);
t50 = cos(t53);
t48 = -t49 * t54 + t50 * t64;
t47 = 0.1e1 / t48 ^ 2;
t70 = -t73 * t90 - t76 * t80;
t85 = t70 * t77 + t73 * t93;
t71 = -t73 * t91 + t76 * t82;
t95 = t71 * t79;
t57 = -t85 * t81 + t95;
t97 = t47 * t57;
t61 = 0.1e1 / t64 ^ 2;
t96 = t54 * t61;
t89 = t79 * t82;
t88 = t80 * t81;
t86 = -t49 * t64 - t50 * t54;
t68 = (t77 * t88 + t89) * t75;
t66 = t73 * t75 * t77 - t70 * t74;
t65 = t79 * t94 + (t77 * t89 + t88) * t75;
t63 = 0.1e1 / t66 ^ 2;
t62 = 0.1e1 / t66;
t60 = 0.1e1 / t64;
t59 = t69 * t92 + t84 * t79;
t58 = t71 * t81 + t85 * t79;
t56 = t69 * t81 + t98 * t79;
t52 = 0.1e1 / (t58 ^ 2 * t63 + 0.1e1);
t51 = 0.1e1 / (t54 ^ 2 * t61 + 0.1e1);
t46 = 0.1e1 / t48;
t45 = 0.1e1 / (t57 ^ 2 * t47 + 0.1e1);
t44 = (-t59 * t60 + t68 * t96) * t51;
t43 = (-t56 * t60 + t65 * t96) * t51;
t1 = [0, t44, t43, 0, 0, 0; 0 ((t70 * t79 + t71 * t92) * t46 - (t86 * t44 - t49 * t59 + t50 * t68) * t97) * t45 (t58 * t46 - (t86 * t43 - t49 * t56 + t50 * t65) * t97) * t45, 0, 0, 0; 0 ((t70 * t81 - t77 * t95) * t62 - t71 * t74 * t58 * t63) * t52, -t57 * t62 * t52, 0, 0, 0;];
Ja_rot  = t1;
