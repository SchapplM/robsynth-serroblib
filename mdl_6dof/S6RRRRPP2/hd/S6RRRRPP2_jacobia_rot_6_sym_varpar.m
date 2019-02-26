% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP2
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
% Datum: 2019-02-26 22:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:25:31
% EndTime: 2019-02-26 22:25:31
% DurationCPUTime: 0.11s
% Computational Cost: add. (290->20), mult. (335->54), div. (80->9), fcn. (494->9), ass. (0->38)
t69 = qJ(2) + qJ(3);
t68 = cos(t69);
t67 = sin(t69);
t71 = sin(qJ(1));
t78 = t71 * t67;
t63 = atan2(t78, t68);
t60 = sin(t63);
t61 = cos(t63);
t52 = t60 * t78 + t61 * t68;
t51 = 0.1e1 / t52 ^ 2;
t73 = cos(qJ(1));
t85 = t51 * t73 ^ 2;
t72 = cos(qJ(4));
t74 = t73 * t72;
t70 = sin(qJ(4));
t77 = t71 * t70;
t59 = t68 * t74 + t77;
t57 = 0.1e1 / t59 ^ 2;
t75 = t73 * t70;
t76 = t71 * t72;
t58 = -t68 * t75 + t76;
t84 = t58 ^ 2 * t57;
t83 = t57 * t58;
t82 = t60 * t68;
t64 = t67 ^ 2;
t81 = t64 / t68 ^ 2;
t80 = t67 * t73;
t62 = 0.1e1 / (t71 ^ 2 * t81 + 0.1e1);
t79 = t71 * t62;
t65 = 0.1e1 / t68;
t56 = 0.1e1 / t59;
t54 = 0.1e1 / (0.1e1 + t84);
t53 = (0.1e1 + t81) * t79;
t50 = 0.1e1 / t52;
t49 = 0.1e1 / (t64 * t85 + 0.1e1);
t48 = (t56 * t70 + t72 * t83) * t54 * t80;
t47 = (-t68 * t50 + (t71 * t82 - t61 * t67 + (t61 * t78 - t82) * t53) * t67 * t51) * t73 * t49;
t1 = [t65 * t62 * t80, t53, t53, 0, 0, 0; (t50 * t78 + (t61 * t64 * t65 * t79 + (-t62 + 0.1e1) * t67 * t60) * t67 * t85) * t49, t47, t47, 0, 0, 0; ((t68 * t77 + t74) * t56 - (-t68 * t76 + t75) * t83) * t54, t48, t48 (-t56 * t59 - t84) * t54, 0, 0;];
Ja_rot  = t1;
