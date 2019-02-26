% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:10:19
% EndTime: 2019-02-26 21:10:19
% DurationCPUTime: 0.13s
% Computational Cost: add. (740->27), mult. (688->69), div. (136->11), fcn. (1034->9), ass. (0->43)
t80 = pkin(10) + qJ(3) + qJ(4);
t76 = sin(t80);
t100 = t76 ^ 2;
t77 = cos(t80);
t86 = cos(qJ(5));
t87 = cos(qJ(1));
t89 = t87 * t86;
t84 = sin(qJ(5));
t85 = sin(qJ(1));
t92 = t85 * t84;
t68 = t77 * t92 + t89;
t94 = t76 * t84;
t65 = atan2(-t68, t94);
t61 = sin(t65);
t62 = cos(t65);
t60 = -t61 * t68 + t62 * t94;
t59 = 0.1e1 / t60 ^ 2;
t90 = t87 * t84;
t91 = t85 * t86;
t71 = t77 * t90 - t91;
t99 = t59 * t71;
t97 = t62 * t68;
t96 = t71 ^ 2 * t59;
t74 = 0.1e1 / t76;
t81 = 0.1e1 / t84;
t95 = t74 * t81;
t93 = t76 * t87;
t72 = t77 * t89 + t92;
t67 = 0.1e1 / t72 ^ 2;
t88 = t87 ^ 2 * t100 * t67;
t82 = 0.1e1 / t84 ^ 2;
t75 = 0.1e1 / t100;
t70 = t77 * t91 - t90;
t66 = 0.1e1 / t72;
t64 = 0.1e1 / (t68 ^ 2 * t75 * t82 + 0.1e1);
t63 = 0.1e1 / (0.1e1 + t88);
t58 = 0.1e1 / t60;
t57 = (t68 * t75 * t77 * t81 + t85) * t64;
t56 = 0.1e1 / (0.1e1 + t96);
t55 = (t68 * t82 * t86 - t70 * t81) * t74 * t64;
t54 = (-t66 * t77 * t87 - t86 * t88) * t63;
t53 = (t57 * t97 * t99 + (-t58 * t93 - (t62 * t77 + (-t57 + t85) * t76 * t61) * t99) * t84) * t56;
t1 = [-t71 * t64 * t95, 0, t57, t57, t55, 0; (-t68 * t58 - (-t61 + (t95 * t97 + t61) * t64) * t96) * t56, 0, t53, t53 (t72 * t58 - (t62 * t76 * t86 - t61 * t70 + (-t61 * t94 - t97) * t55) * t99) * t56, 0; (-t67 * t70 * t87 + t66 * t85) * t76 * t63, 0, t54, t54, -t71 * t67 * t63 * t93, 0;];
Ja_rot  = t1;
