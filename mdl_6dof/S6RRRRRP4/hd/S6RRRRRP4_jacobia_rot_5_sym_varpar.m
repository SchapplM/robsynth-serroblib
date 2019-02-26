% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRP4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:41:20
% EndTime: 2019-02-26 22:41:20
% DurationCPUTime: 0.11s
% Computational Cost: add. (453->22), mult. (359->54), div. (84->9), fcn. (524->9), ass. (0->40)
t80 = qJ(2) + qJ(3);
t78 = cos(t80);
t76 = sin(t80);
t81 = sin(qJ(1));
t87 = t81 * t76;
t70 = atan2(-t87, -t78);
t68 = sin(t70);
t69 = cos(t70);
t61 = -t68 * t87 - t69 * t78;
t60 = 0.1e1 / t61 ^ 2;
t82 = cos(qJ(1));
t94 = t60 * t82 ^ 2;
t79 = qJ(4) + qJ(5);
t77 = cos(t79);
t84 = t82 * t77;
t75 = sin(t79);
t88 = t81 * t75;
t67 = t78 * t84 + t88;
t65 = 0.1e1 / t67 ^ 2;
t85 = t82 * t75;
t86 = t81 * t77;
t66 = t78 * t85 - t86;
t93 = t65 * t66;
t92 = t68 * t78;
t72 = t76 ^ 2;
t91 = t72 / t78 ^ 2;
t90 = t76 * t82;
t71 = 0.1e1 / (t81 ^ 2 * t91 + 0.1e1);
t89 = t81 * t71;
t83 = t66 ^ 2 * t65 + 0.1e1;
t73 = 0.1e1 / t78;
t64 = 0.1e1 / t67;
t63 = (0.1e1 + t91) * t89;
t62 = 0.1e1 / t83;
t59 = 0.1e1 / t61;
t58 = 0.1e1 / (t72 * t94 + 0.1e1);
t57 = t83 * t62;
t56 = (-t64 * t75 + t77 * t93) * t62 * t90;
t55 = (t78 * t59 - (-t81 * t92 + t69 * t76 + (-t69 * t87 + t92) * t63) * t76 * t60) * t82 * t58;
t1 = [t73 * t71 * t90, t63, t63, 0, 0, 0; (-t59 * t87 - (-t69 * t72 * t73 * t89 + (t71 - 0.1e1) * t76 * t68) * t76 * t94) * t58, t55, t55, 0, 0, 0; ((-t78 * t88 - t84) * t64 - (-t78 * t86 + t85) * t93) * t62, t56, t56, t57, t57, 0;];
Ja_rot  = t1;
