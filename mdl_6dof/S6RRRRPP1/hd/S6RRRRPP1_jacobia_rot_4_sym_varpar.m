% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP1_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_jacobia_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:58
% EndTime: 2019-02-26 22:24:58
% DurationCPUTime: 0.11s
% Computational Cost: add. (358->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
t68 = qJ(2) + qJ(3);
t67 = cos(t68);
t66 = sin(t68);
t70 = sin(qJ(1));
t78 = t70 * t66;
t61 = atan2(-t78, -t67);
t59 = sin(t61);
t60 = cos(t61);
t52 = -t59 * t78 - t60 * t67;
t51 = 0.1e1 / t52 ^ 2;
t72 = cos(qJ(1));
t84 = t51 * t72 ^ 2;
t71 = cos(qJ(4));
t74 = t72 * t71;
t69 = sin(qJ(4));
t77 = t70 * t69;
t58 = t67 * t74 + t77;
t56 = 0.1e1 / t58 ^ 2;
t75 = t72 * t69;
t76 = t70 * t71;
t57 = t67 * t75 - t76;
t83 = t56 * t57;
t82 = t59 * t67;
t63 = t66 ^ 2;
t81 = t63 / t67 ^ 2;
t80 = t66 * t72;
t62 = 0.1e1 / (t70 ^ 2 * t81 + 0.1e1);
t79 = t70 * t62;
t73 = t57 ^ 2 * t56 + 0.1e1;
t64 = 0.1e1 / t67;
t55 = 0.1e1 / t58;
t54 = 0.1e1 / t73;
t53 = (0.1e1 + t81) * t79;
t50 = 0.1e1 / t52;
t49 = 0.1e1 / (t63 * t84 + 0.1e1);
t48 = (-t55 * t69 + t71 * t83) * t54 * t80;
t47 = (t67 * t50 - (-t70 * t82 + t60 * t66 + (-t60 * t78 + t82) * t53) * t66 * t51) * t72 * t49;
t1 = [t64 * t62 * t80, t53, t53, 0, 0, 0; (-t50 * t78 - (-t60 * t63 * t64 * t79 + (t62 - 0.1e1) * t66 * t59) * t66 * t84) * t49, t47, t47, 0, 0, 0; ((-t67 * t77 - t74) * t55 - (-t67 * t76 + t75) * t83) * t54, t48, t48, t73 * t54, 0, 0;];
Ja_rot  = t1;
