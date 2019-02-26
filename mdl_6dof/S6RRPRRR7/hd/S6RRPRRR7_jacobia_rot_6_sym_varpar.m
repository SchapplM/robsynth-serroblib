% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:57:35
% EndTime: 2019-02-26 21:57:35
% DurationCPUTime: 0.20s
% Computational Cost: add. (402->24), mult. (924->55), div. (90->9), fcn. (1345->11), ass. (0->39)
t102 = sin(qJ(2));
t103 = cos(qJ(4));
t84 = sin(qJ(4));
t86 = cos(qJ(2));
t75 = t102 * t103 - t86 * t84;
t85 = sin(qJ(1));
t67 = t75 * t85;
t74 = t102 * t84 + t86 * t103;
t73 = 0.1e1 / t74 ^ 2;
t61 = 0.1e1 / (t67 ^ 2 * t73 + 0.1e1);
t69 = t74 * t85;
t72 = 0.1e1 / t74;
t106 = (-t67 * t73 * t75 - t69 * t72) * t61;
t64 = atan2(t67, t74);
t59 = sin(t64);
t60 = cos(t64);
t57 = t59 * t67 + t60 * t74;
t56 = 0.1e1 / t57 ^ 2;
t104 = cos(qJ(1));
t71 = t75 * t104;
t96 = t71 ^ 2 * t56;
t54 = 0.1e1 / (0.1e1 + t96);
t55 = 0.1e1 / t57;
t70 = t74 * t104;
t99 = t60 * t67;
t109 = ((t106 * (-t59 * t74 + t99) - t59 * t69 + t60 * t75) * t56 * t71 + t70 * t55) * t54;
t83 = qJ(5) + qJ(6);
t81 = sin(t83);
t82 = cos(t83);
t66 = t70 * t82 - t85 * t81;
t63 = 0.1e1 / t66 ^ 2;
t65 = t70 * t81 + t85 * t82;
t94 = t65 ^ 2 * t63 + 0.1e1;
t58 = 0.1e1 / t94;
t62 = 0.1e1 / t66;
t98 = t63 * t65;
t105 = (-t81 * t62 + t82 * t98) * t58 * t71;
t51 = t94 * t58;
t1 = [t71 * t72 * t61, -t106, 0, t106, 0, 0; (t67 * t55 - (-t59 + (-t72 * t99 + t59) * t61) * t96) * t54, -t109, 0, t109, 0, 0; ((t104 * t82 - t69 * t81) * t62 - (-t104 * t81 - t69 * t82) * t98) * t58, t105, 0, -t105, t51, t51;];
Ja_rot  = t1;
