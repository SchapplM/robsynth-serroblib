% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRPRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRPRR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:48
% EndTime: 2019-02-26 19:39:48
% DurationCPUTime: 0.19s
% Computational Cost: add. (582->31), mult. (1659->71), div. (35->9), fcn. (2273->17), ass. (0->48)
t83 = sin(pkin(11));
t85 = sin(pkin(6));
t100 = t83 * t85;
t90 = cos(pkin(6));
t99 = t83 * t90;
t88 = cos(pkin(11));
t98 = t85 * t88;
t97 = t88 * t90;
t84 = sin(pkin(7));
t81 = sin(pkin(13));
t86 = cos(pkin(13));
t92 = sin(qJ(3));
t94 = cos(qJ(3));
t95 = t94 * t81 + t92 * t86;
t72 = t95 * t84;
t89 = cos(pkin(7));
t74 = t95 * t89;
t82 = sin(pkin(12));
t87 = cos(pkin(12));
t77 = -t88 * t82 - t87 * t99;
t78 = -t82 * t99 + t88 * t87;
t79 = t92 * t81 - t94 * t86;
t65 = t72 * t100 + t77 * t74 - t78 * t79;
t69 = t89 * t100 - t77 * t84;
t91 = sin(qJ(5));
t93 = cos(qJ(5));
t60 = t65 * t93 + t69 * t91;
t58 = 0.1e1 / t60 ^ 2;
t59 = t65 * t91 - t69 * t93;
t96 = t59 ^ 2 * t58 + 0.1e1;
t76 = t82 * t97 + t83 * t87;
t75 = -t83 * t82 + t87 * t97;
t73 = t79 * t89;
t71 = t79 * t84;
t68 = t90 * t71 + (t73 * t87 + t82 * t95) * t85;
t67 = t90 * t72 + (t74 * t87 - t79 * t82) * t85;
t66 = 0.1e1 / t68 ^ 2;
t63 = -t71 * t100 - t77 * t73 - t78 * t95;
t62 = t71 * t98 - t75 * t73 - t76 * t95;
t61 = t72 * t98 - t75 * t74 + t76 * t79;
t57 = atan2(t62, t68);
t55 = cos(t57);
t54 = sin(t57);
t53 = 0.1e1 / t96;
t52 = t54 * t62 + t55 * t68;
t51 = 0.1e1 / t52 ^ 2;
t49 = (t61 / t68 - t67 * t62 * t66) / (t62 ^ 2 * t66 + 0.1e1);
t1 = [0, 0, t49, 0, 0, 0; 0, 0 (t65 / t52 + (t54 * t61 + t55 * t67 + (-t54 * t68 + t55 * t62) * t49) * t63 * t51) / (t63 ^ 2 * t51 + 0.1e1) 0, 0, 0; 0, 0 (t91 / t60 - t93 * t59 * t58) * t63 * t53, 0, t96 * t53, 0;];
Ja_rot  = t1;
