% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:31:16
% EndTime: 2019-02-26 22:31:16
% DurationCPUTime: 0.12s
% Computational Cost: add. (717->21), mult. (413->53), div. (101->9), fcn. (611->9), ass. (0->37)
t68 = qJ(2) + qJ(3) + qJ(4);
t67 = cos(t68);
t66 = sin(t68);
t71 = sin(qJ(1));
t77 = t71 * t66;
t61 = atan2(-t77, -t67);
t57 = sin(t61);
t58 = cos(t61);
t52 = -t57 * t77 - t58 * t67;
t51 = 0.1e1 / t52 ^ 2;
t72 = cos(qJ(1));
t83 = t51 * t72 ^ 2;
t70 = cos(pkin(11));
t73 = t72 * t70;
t69 = sin(pkin(11));
t76 = t71 * t69;
t60 = t67 * t73 + t76;
t56 = 0.1e1 / t60 ^ 2;
t74 = t72 * t69;
t75 = t71 * t70;
t59 = t67 * t74 - t75;
t82 = t56 * t59;
t81 = t57 * t67;
t63 = t66 ^ 2;
t80 = t63 / t67 ^ 2;
t79 = t66 * t72;
t62 = 0.1e1 / (t71 ^ 2 * t80 + 0.1e1);
t78 = t71 * t62;
t64 = 0.1e1 / t67;
t55 = 0.1e1 / t60;
t54 = 0.1e1 / (t59 ^ 2 * t56 + 0.1e1);
t53 = (0.1e1 + t80) * t78;
t50 = 0.1e1 / t52;
t49 = 0.1e1 / (t63 * t83 + 0.1e1);
t48 = (-t55 * t69 + t70 * t82) * t54 * t79;
t47 = (t67 * t50 - (-t71 * t81 + t58 * t66 + (-t58 * t77 + t81) * t53) * t66 * t51) * t72 * t49;
t1 = [t64 * t62 * t79, t53, t53, t53, 0, 0; (-t50 * t77 - (-t58 * t63 * t64 * t78 + (t62 - 0.1e1) * t66 * t57) * t66 * t83) * t49, t47, t47, t47, 0, 0; ((-t67 * t76 - t73) * t55 - (-t67 * t75 + t74) * t82) * t54, t48, t48, t48, 0, 0;];
Ja_rot  = t1;
