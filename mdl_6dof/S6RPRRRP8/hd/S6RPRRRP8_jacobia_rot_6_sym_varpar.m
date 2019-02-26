% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP8
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
% Datum: 2019-02-26 21:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:12:04
% EndTime: 2019-02-26 21:12:04
% DurationCPUTime: 0.14s
% Computational Cost: add. (478->25), mult. (688->69), div. (136->11), fcn. (1034->9), ass. (0->43)
t74 = qJ(3) + qJ(4);
t70 = cos(t74);
t91 = t70 ^ 2;
t69 = sin(t74);
t75 = sin(qJ(5));
t78 = cos(qJ(1));
t81 = t78 * t75;
t76 = sin(qJ(1));
t77 = cos(qJ(5));
t82 = t76 * t77;
t63 = t69 * t81 + t82;
t85 = t70 * t75;
t58 = atan2(t63, t85);
t54 = sin(t58);
t55 = cos(t58);
t53 = t54 * t63 + t55 * t85;
t52 = 0.1e1 / t53 ^ 2;
t80 = t78 * t77;
t83 = t76 * t75;
t61 = t69 * t83 - t80;
t90 = t52 * t61;
t88 = t55 * t63;
t87 = t61 ^ 2 * t52;
t67 = 0.1e1 / t70;
t71 = 0.1e1 / t75;
t86 = t67 * t71;
t84 = t70 * t76;
t62 = t69 * t82 + t81;
t60 = 0.1e1 / t62 ^ 2;
t79 = t76 ^ 2 * t91 * t60;
t72 = 0.1e1 / t75 ^ 2;
t68 = 0.1e1 / t91;
t64 = t69 * t80 - t83;
t59 = 0.1e1 / t62;
t57 = 0.1e1 / (t63 ^ 2 * t68 * t72 + 0.1e1);
t56 = 0.1e1 / (0.1e1 + t79);
t51 = 0.1e1 / t53;
t50 = (t63 * t68 * t69 * t71 + t78) * t57;
t49 = 0.1e1 / (0.1e1 + t87);
t48 = (-t63 * t72 * t77 + t64 * t71) * t67 * t57;
t47 = (-t59 * t69 * t76 - t77 * t79) * t56;
t46 = (-t50 * t88 * t90 + (t51 * t84 - (-t55 * t69 + (-t50 + t78) * t54 * t70) * t90) * t75) * t49;
t1 = [-t61 * t57 * t86, 0, t50, t50, t48, 0; (t63 * t51 - (-t54 + (-t86 * t88 + t54) * t57) * t87) * t49, 0, t46, t46 (t62 * t51 - (t55 * t70 * t77 + t54 * t64 + (-t54 * t85 + t88) * t48) * t90) * t49, 0; (-t60 * t64 * t76 + t59 * t78) * t70 * t56, 0, t47, t47, t61 * t60 * t56 * t84, 0;];
Ja_rot  = t1;
