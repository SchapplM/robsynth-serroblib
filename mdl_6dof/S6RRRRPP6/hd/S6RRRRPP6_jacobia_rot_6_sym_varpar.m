% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP6
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
% Datum: 2019-02-26 22:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:28:14
% EndTime: 2019-02-26 22:28:14
% DurationCPUTime: 0.14s
% Computational Cost: add. (598->28), mult. (659->69), div. (149->11), fcn. (1026->9), ass. (0->41)
t71 = qJ(3) + qJ(4);
t66 = cos(t71);
t65 = sin(t71);
t75 = cos(qJ(1));
t77 = t75 * t65;
t73 = sin(qJ(1));
t74 = cos(qJ(2));
t78 = t73 * t74;
t58 = t66 * t78 - t77;
t72 = sin(qJ(2));
t79 = t72 * t66;
t55 = atan2(-t58, t79);
t52 = sin(t55);
t53 = cos(t55);
t50 = -t52 * t58 + t53 * t79;
t49 = 0.1e1 / t50 ^ 2;
t76 = t75 * t66;
t61 = t73 * t65 + t74 * t76;
t86 = t49 * t61;
t84 = t53 * t58;
t60 = t73 * t66 - t74 * t77;
t68 = 0.1e1 / t72 ^ 2;
t70 = 0.1e1 / t75 ^ 2;
t56 = 0.1e1 / (t60 ^ 2 * t70 * t68 + 0.1e1);
t67 = 0.1e1 / t72;
t83 = t56 * t67;
t82 = t61 ^ 2 * t49;
t63 = 0.1e1 / t66;
t81 = t63 * t67;
t80 = t68 * t74;
t69 = 0.1e1 / t75;
t64 = 0.1e1 / t66 ^ 2;
t57 = t65 * t78 + t76;
t54 = 0.1e1 / (t58 ^ 2 * t68 * t64 + 0.1e1);
t51 = t61 * t69 * t83;
t48 = 0.1e1 / t50;
t47 = (t58 * t63 * t80 + t73) * t54;
t46 = 0.1e1 / (0.1e1 + t82);
t45 = (-t58 * t64 * t65 + t57 * t63) * t67 * t54;
t44 = (t60 * t48 - (-t53 * t72 * t65 + t52 * t57 + (-t52 * t79 - t84) * t45) * t86) * t46;
t1 = [-t61 * t54 * t81, t47, t45, t45, 0, 0; (-t58 * t48 - (-t52 + (t81 * t84 + t52) * t54) * t82) * t46 (t47 * t84 * t86 + (-t75 * t72 * t48 - (t53 * t74 + (-t47 + t73) * t52 * t72) * t86) * t66) * t46, t44, t44, 0, 0; (t60 * t70 * t73 + t57 * t69) * t83 (-t60 * t69 * t80 + t65) * t56, -t51, -t51, 0, 0;];
Ja_rot  = t1;
