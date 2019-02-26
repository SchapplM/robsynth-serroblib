% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:46
% EndTime: 2019-02-26 21:55:46
% DurationCPUTime: 0.19s
% Computational Cost: add. (467->27), mult. (1134->69), div. (60->9), fcn. (1602->13), ass. (0->47)
t89 = sin(qJ(1));
t91 = cos(qJ(1));
t85 = sin(pkin(12));
t87 = cos(pkin(12));
t90 = cos(qJ(2));
t98 = cos(pkin(6));
t95 = t90 * t98;
t88 = sin(qJ(2));
t96 = t88 * t98;
t92 = -t85 * t96 + t87 * t95;
t93 = t90 * t85 + t88 * t87;
t68 = -t89 * t93 + t91 * t92;
t79 = t88 * t85 - t90 * t87;
t86 = sin(pkin(6));
t76 = t79 * t86;
t62 = atan2(t68, t76);
t60 = cos(t62);
t103 = t60 * t68;
t100 = t86 * t89;
t78 = t85 * t95 + t87 * t96;
t72 = -t89 * t78 - t91 * t79;
t84 = qJ(4) + qJ(5);
t82 = sin(t84);
t83 = cos(t84);
t66 = t82 * t100 + t72 * t83;
t64 = 0.1e1 / t66 ^ 2;
t65 = -t83 * t100 + t72 * t82;
t102 = t64 * t65;
t59 = sin(t62);
t57 = t59 * t68 + t60 * t76;
t56 = 0.1e1 / t57 ^ 2;
t70 = -t89 * t92 - t91 * t93;
t101 = t70 ^ 2 * t56;
t99 = t86 * t91;
t97 = t65 ^ 2 * t64 + 0.1e1;
t94 = -t91 * t78 + t89 * t79;
t77 = t93 * t86;
t75 = 0.1e1 / t76 ^ 2;
t74 = 0.1e1 / t76;
t63 = 0.1e1 / t66;
t61 = 0.1e1 / (t68 ^ 2 * t75 + 0.1e1);
t58 = 0.1e1 / t97;
t55 = 0.1e1 / t57;
t54 = 0.1e1 / (0.1e1 + t101);
t53 = (-t68 * t75 * t77 + t74 * t94) * t61;
t52 = t97 * t58;
t1 = [t70 * t74 * t61, t53, 0, 0, 0, 0; (t68 * t55 + (t59 + (t74 * t103 - t59) * t61) * t101) * t54 (t72 * t55 + (t59 * t94 + t60 * t77 + (-t59 * t76 + t103) * t53) * t70 * t56) * t54, 0, 0, 0, 0; ((t82 * t94 - t83 * t99) * t63 - (t82 * t99 + t83 * t94) * t102) * t58 (-t83 * t102 + t82 * t63) * t70 * t58, 0, t52, t52, 0;];
Ja_rot  = t1;
