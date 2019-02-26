% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:26:19
% EndTime: 2019-02-26 22:26:20
% DurationCPUTime: 0.13s
% Computational Cost: add. (469->27), mult. (656->69), div. (146->11), fcn. (1013->9), ass. (0->42)
t82 = qJ(2) + qJ(3);
t77 = cos(t82);
t85 = cos(qJ(4));
t86 = cos(qJ(1));
t87 = t86 * t85;
t83 = sin(qJ(4));
t84 = sin(qJ(1));
t90 = t84 * t83;
t67 = t77 * t90 + t87;
t76 = sin(t82);
t91 = t76 * t83;
t66 = atan2(-t67, t91);
t62 = sin(t66);
t63 = cos(t66);
t61 = -t62 * t67 + t63 * t91;
t60 = 0.1e1 / t61 ^ 2;
t88 = t86 * t83;
t89 = t84 * t85;
t70 = t77 * t88 - t89;
t98 = t60 * t70;
t96 = t63 * t67;
t71 = t77 * t87 + t90;
t75 = 0.1e1 / t76 ^ 2;
t81 = 0.1e1 / t86 ^ 2;
t65 = 0.1e1 / (t71 ^ 2 * t81 * t75 + 0.1e1);
t74 = 0.1e1 / t76;
t95 = t65 * t74;
t94 = t70 ^ 2 * t60;
t78 = 0.1e1 / t83;
t93 = t74 * t78;
t92 = t75 * t77;
t80 = 0.1e1 / t86;
t79 = 0.1e1 / t83 ^ 2;
t69 = t77 * t89 - t88;
t64 = 0.1e1 / (t67 ^ 2 * t75 * t79 + 0.1e1);
t59 = 0.1e1 / t61;
t58 = (-t71 * t80 * t92 - t85) * t65;
t57 = (t67 * t78 * t92 + t84) * t64;
t56 = 0.1e1 / (0.1e1 + t94);
t55 = (t67 * t79 * t85 - t69 * t78) * t74 * t64;
t54 = (t57 * t96 * t98 + (-t86 * t76 * t59 - (t63 * t77 + (-t57 + t84) * t76 * t62) * t98) * t83) * t56;
t1 = [-t70 * t64 * t93, t57, t57, t55, 0, 0; (-t67 * t59 - (-t62 + (t93 * t96 + t62) * t64) * t94) * t56, t54, t54 (t71 * t59 - (t63 * t76 * t85 - t62 * t69 + (-t62 * t91 - t96) * t55) * t98) * t56, 0, 0; (t71 * t81 * t84 - t69 * t80) * t95, t58, t58, -t70 * t80 * t95, 0, 0;];
Ja_rot  = t1;
