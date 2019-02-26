% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP10_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:04
% EndTime: 2019-02-26 21:13:04
% DurationCPUTime: 0.14s
% Computational Cost: add. (612->25), mult. (689->67), div. (139->11), fcn. (1046->9), ass. (0->43)
t76 = cos(qJ(3));
t91 = t76 ^ 2;
t74 = sin(qJ(3));
t73 = qJ(4) + qJ(5);
t67 = sin(t73);
t77 = cos(qJ(1));
t81 = t77 * t67;
t68 = cos(t73);
t75 = sin(qJ(1));
t84 = t75 * t68;
t62 = t74 * t81 + t84;
t82 = t76 * t67;
t56 = atan2(t62, t82);
t53 = sin(t56);
t54 = cos(t56);
t52 = t53 * t62 + t54 * t82;
t51 = 0.1e1 / t52 ^ 2;
t80 = t77 * t68;
t85 = t75 * t67;
t60 = t74 * t85 - t80;
t90 = t51 * t60;
t88 = t54 * t62;
t87 = t60 ^ 2 * t51;
t65 = 0.1e1 / t67;
t71 = 0.1e1 / t76;
t86 = t65 * t71;
t83 = t75 * t76;
t61 = t74 * t84 + t81;
t59 = 0.1e1 / t61 ^ 2;
t79 = t75 ^ 2 * t91 * t59;
t57 = 0.1e1 / (0.1e1 + t79);
t78 = t60 * t59 * t57 * t83;
t72 = 0.1e1 / t91;
t66 = 0.1e1 / t67 ^ 2;
t63 = t74 * t80 - t85;
t58 = 0.1e1 / t61;
t55 = 0.1e1 / (t62 ^ 2 * t66 * t72 + 0.1e1);
t50 = 0.1e1 / t52;
t49 = (t62 * t65 * t72 * t74 + t77) * t55;
t48 = 0.1e1 / (0.1e1 + t87);
t47 = (-t62 * t66 * t68 + t63 * t65) * t71 * t55;
t46 = (t61 * t50 - (t54 * t76 * t68 + t53 * t63 + (-t53 * t82 + t88) * t47) * t90) * t48;
t1 = [-t60 * t55 * t86, 0, t49, t47, t47, 0; (t62 * t50 - (-t53 + (-t86 * t88 + t53) * t55) * t87) * t48, 0 (-t49 * t88 * t90 + (t50 * t83 - (-t54 * t74 + (-t49 + t77) * t76 * t53) * t90) * t67) * t48, t46, t46, 0; (-t59 * t63 * t75 + t58 * t77) * t76 * t57, 0 (-t58 * t74 * t75 - t68 * t79) * t57, t78, t78, 0;];
Ja_rot  = t1;
