% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:37:04
% EndTime: 2019-02-26 20:37:04
% DurationCPUTime: 0.11s
% Computational Cost: add. (288->19), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
t73 = qJ(4) + qJ(5);
t71 = sin(t73);
t72 = cos(t73);
t75 = sin(qJ(1));
t83 = t75 * t72;
t67 = atan2(t83, t71);
t64 = sin(t67);
t65 = cos(t67);
t57 = t64 * t83 + t65 * t71;
t56 = 0.1e1 / t57 ^ 2;
t77 = cos(qJ(1));
t89 = t56 * t77 ^ 2;
t76 = cos(qJ(6));
t79 = t77 * t76;
t74 = sin(qJ(6));
t82 = t75 * t74;
t63 = t71 * t79 - t82;
t61 = 0.1e1 / t63 ^ 2;
t80 = t77 * t74;
t81 = t75 * t76;
t62 = t71 * t80 + t81;
t88 = t61 * t62;
t87 = t64 * t71;
t70 = t72 ^ 2;
t86 = 0.1e1 / t71 ^ 2 * t70;
t85 = t72 * t77;
t66 = 0.1e1 / (t75 ^ 2 * t86 + 0.1e1);
t84 = t75 * t66;
t78 = t62 ^ 2 * t61 + 0.1e1;
t68 = 0.1e1 / t71;
t60 = 0.1e1 / t63;
t59 = 0.1e1 / t78;
t58 = (-0.1e1 - t86) * t84;
t55 = 0.1e1 / t57;
t54 = 0.1e1 / (t70 * t89 + 0.1e1);
t53 = (t60 * t74 - t76 * t88) * t59 * t85;
t52 = (t71 * t55 + (-t75 * t87 + t65 * t72 + (t65 * t83 - t87) * t58) * t72 * t56) * t77 * t54;
t1 = [t68 * t66 * t85, 0, 0, t58, t58, 0; (t55 * t83 + (t65 * t68 * t70 * t84 + (-t66 + 0.1e1) * t72 * t64) * t72 * t89) * t54, 0, 0, t52, t52, 0; ((-t71 * t82 + t79) * t60 - (-t71 * t81 - t80) * t88) * t59, 0, 0, t53, t53, t78 * t59;];
Ja_rot  = t1;
