% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:07
% EndTime: 2019-02-26 20:11:07
% DurationCPUTime: 0.19s
% Computational Cost: add. (903->31), mult. (1304->78), div. (90->9), fcn. (1864->13), ass. (0->53)
t89 = sin(pkin(6));
t91 = cos(pkin(11));
t100 = t89 * t91;
t88 = sin(pkin(11));
t94 = cos(qJ(2));
t92 = cos(pkin(6));
t93 = sin(qJ(2));
t97 = t92 * t93;
t80 = t88 * t94 + t91 * t97;
t86 = qJ(3) + qJ(4);
t84 = sin(t86);
t85 = cos(t86);
t70 = t85 * t100 + t80 * t84;
t99 = t89 * t93;
t77 = t84 * t99 - t92 * t85;
t69 = atan2(-t70, t77);
t66 = sin(t69);
t67 = cos(t69);
t60 = -t66 * t70 + t67 * t77;
t59 = 0.1e1 / t60 ^ 2;
t101 = t88 * t89;
t82 = -t88 * t97 + t91 * t94;
t73 = -t85 * t101 + t82 * t84;
t106 = t59 * t73;
t96 = t92 * t94;
t81 = t88 * t96 + t91 * t93;
t87 = sin(pkin(12));
t103 = t81 * t87;
t74 = t84 * t101 + t82 * t85;
t90 = cos(pkin(12));
t65 = t74 * t90 + t103;
t63 = 0.1e1 / t65 ^ 2;
t102 = t81 * t90;
t64 = t74 * t87 - t102;
t105 = t63 * t64;
t76 = 0.1e1 / t77 ^ 2;
t104 = t70 * t76;
t98 = t89 * t94;
t95 = -t66 * t77 - t67 * t70;
t79 = -t88 * t93 + t91 * t96;
t78 = t92 * t84 + t85 * t99;
t75 = 0.1e1 / t77;
t72 = -t84 * t100 + t80 * t85;
t68 = 0.1e1 / (t70 ^ 2 * t76 + 0.1e1);
t62 = 0.1e1 / t65;
t61 = 0.1e1 / (t64 ^ 2 * t63 + 0.1e1);
t58 = 0.1e1 / t60;
t57 = 0.1e1 / (t73 ^ 2 * t59 + 0.1e1);
t56 = (t98 * t104 - t75 * t79) * t84 * t68;
t55 = (t78 * t104 - t72 * t75) * t68;
t54 = (t90 * t105 - t62 * t87) * t73 * t61;
t53 = (t74 * t58 - (t95 * t55 - t66 * t72 + t67 * t78) * t106) * t57;
t1 = [0, t56, t55, t55, 0, 0; 0 (-t81 * t84 * t58 - ((-t66 * t79 + t67 * t98) * t84 + t95 * t56) * t106) * t57, t53, t53, 0, 0; 0 ((-t85 * t103 - t82 * t90) * t62 - (-t85 * t102 + t82 * t87) * t105) * t61, t54, t54, 0, 0;];
Ja_rot  = t1;
