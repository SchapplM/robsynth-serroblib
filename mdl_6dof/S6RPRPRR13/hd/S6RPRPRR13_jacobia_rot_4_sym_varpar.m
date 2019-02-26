% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR13_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobia_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:38
% EndTime: 2019-02-26 20:55:38
% DurationCPUTime: 0.17s
% Computational Cost: add. (423->35), mult. (1287->76), div. (47->9), fcn. (1778->13), ass. (0->50)
t75 = sin(pkin(7));
t78 = cos(pkin(7));
t79 = cos(pkin(6));
t74 = sin(pkin(12));
t83 = cos(qJ(1));
t90 = t83 * t74;
t77 = cos(pkin(12));
t81 = sin(qJ(1));
t91 = t81 * t77;
t86 = t79 * t91 + t90;
t76 = sin(pkin(6));
t96 = t76 * t81;
t101 = -t75 * t96 + t86 * t78;
t89 = t83 * t77;
t92 = t81 * t74;
t69 = -t79 * t89 + t92;
t82 = cos(qJ(3));
t95 = t76 * t83;
t88 = t75 * t95;
t93 = t78 * t82;
t70 = t79 * t90 + t91;
t80 = sin(qJ(3));
t98 = t70 * t80;
t55 = t69 * t93 + t82 * t88 + t98;
t97 = t75 * t79;
t62 = -t82 * t97 + (t74 * t80 - t77 * t93) * t76;
t53 = atan2(-t55, t62);
t51 = cos(t53);
t100 = t51 * t55;
t50 = sin(t53);
t49 = -t50 * t55 + t51 * t62;
t48 = 0.1e1 / t49 ^ 2;
t71 = -t79 * t92 + t89;
t58 = t101 * t82 + t71 * t80;
t99 = t58 ^ 2 * t48;
t94 = t78 * t80;
t85 = t69 * t94 - t70 * t82 + t80 * t88;
t66 = t86 * t75 + t78 * t96;
t65 = 0.1e1 / t66 ^ 2;
t64 = 0.1e1 / t66;
t63 = t80 * t97 + (t74 * t82 + t77 * t94) * t76;
t61 = 0.1e1 / t62 ^ 2;
t60 = 0.1e1 / t62;
t59 = -t101 * t80 + t71 * t82;
t54 = 0.1e1 / (t59 ^ 2 * t65 + 0.1e1);
t52 = 0.1e1 / (t55 ^ 2 * t61 + 0.1e1);
t47 = 0.1e1 / t49;
t46 = 0.1e1 / (0.1e1 + t99);
t45 = (t55 * t61 * t63 + t60 * t85) * t52;
t1 = [-t58 * t60 * t52, 0, t45, 0, 0, 0; ((-t98 + (-t69 * t78 - t88) * t82) * t47 - (-t50 + (t60 * t100 + t50) * t52) * t99) * t46, 0 (t59 * t47 - (t50 * t85 + t51 * t63 + (-t50 * t62 - t100) * t45) * t58 * t48) * t46, 0, 0, 0; (t85 * t64 - (-t69 * t75 + t78 * t95) * t59 * t65) * t54, 0, -t58 * t64 * t54, 0, 0, 0;];
Ja_rot  = t1;
