% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_rot = S6RRRRPP6_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:28:23
% EndTime: 2019-02-26 22:28:24
% DurationCPUTime: 0.13s
% Computational Cost: add. (598->28), mult. (659->69), div. (149->11), fcn. (1026->9), ass. (0->41)
t73 = qJ(3) + qJ(4);
t67 = sin(t73);
t68 = cos(t73);
t77 = cos(qJ(1));
t78 = t77 * t68;
t75 = sin(qJ(1));
t76 = cos(qJ(2));
t80 = t75 * t76;
t58 = t67 * t80 + t78;
t74 = sin(qJ(2));
t81 = t74 * t67;
t56 = atan2(-t58, t81);
t53 = sin(t56);
t54 = cos(t56);
t51 = -t53 * t58 + t54 * t81;
t50 = 0.1e1 / t51 ^ 2;
t79 = t77 * t67;
t61 = -t75 * t68 + t76 * t79;
t88 = t50 * t61;
t86 = t54 * t58;
t62 = t75 * t67 + t76 * t78;
t70 = 0.1e1 / t74 ^ 2;
t72 = 0.1e1 / t77 ^ 2;
t57 = 0.1e1 / (t62 ^ 2 * t72 * t70 + 0.1e1);
t69 = 0.1e1 / t74;
t85 = t57 * t69;
t84 = t61 ^ 2 * t50;
t65 = 0.1e1 / t67;
t83 = t65 * t69;
t82 = t70 * t76;
t71 = 0.1e1 / t77;
t66 = 0.1e1 / t67 ^ 2;
t60 = t68 * t80 - t79;
t55 = 0.1e1 / (t58 ^ 2 * t70 * t66 + 0.1e1);
t52 = t61 * t71 * t85;
t49 = 0.1e1 / t51;
t48 = (t58 * t65 * t82 + t75) * t55;
t47 = 0.1e1 / (0.1e1 + t84);
t46 = (t58 * t66 * t68 - t60 * t65) * t69 * t55;
t45 = (t62 * t49 - (t54 * t74 * t68 - t53 * t60 + (-t53 * t81 - t86) * t46) * t88) * t47;
t1 = [-t61 * t55 * t83, t48, t46, t46, 0, 0; (-t58 * t49 - (-t53 + (t83 * t86 + t53) * t55) * t84) * t47 (t48 * t86 * t88 + (-t77 * t74 * t49 - (t54 * t76 + (-t48 + t75) * t53 * t74) * t88) * t67) * t47, t45, t45, 0, 0; (t62 * t72 * t75 - t60 * t71) * t85 (-t62 * t71 * t82 - t68) * t57, -t52, -t52, 0, 0;];
Ja_rot  = t1;
