% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:09:48
% EndTime: 2019-02-26 21:09:48
% DurationCPUTime: 0.11s
% Computational Cost: add. (554->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
t67 = pkin(10) + qJ(3) + qJ(4);
t66 = cos(t67);
t65 = sin(t67);
t69 = sin(qJ(1));
t77 = t69 * t65;
t58 = atan2(-t77, -t66);
t54 = sin(t58);
t55 = cos(t58);
t51 = -t54 * t77 - t55 * t66;
t50 = 0.1e1 / t51 ^ 2;
t71 = cos(qJ(1));
t83 = t50 * t71 ^ 2;
t82 = t54 * t66;
t70 = cos(qJ(5));
t73 = t71 * t70;
t68 = sin(qJ(5));
t76 = t69 * t68;
t60 = t66 * t73 + t76;
t57 = 0.1e1 / t60 ^ 2;
t74 = t71 * t68;
t75 = t69 * t70;
t59 = t66 * t74 - t75;
t81 = t57 * t59;
t62 = t65 ^ 2;
t80 = t62 / t66 ^ 2;
t79 = t65 * t71;
t61 = 0.1e1 / (t69 ^ 2 * t80 + 0.1e1);
t78 = t69 * t61;
t72 = t59 ^ 2 * t57 + 0.1e1;
t63 = 0.1e1 / t66;
t56 = 0.1e1 / t60;
t53 = 0.1e1 / t72;
t52 = (0.1e1 + t80) * t78;
t49 = 0.1e1 / t51;
t48 = 0.1e1 / (t62 * t83 + 0.1e1);
t47 = (-t56 * t68 + t70 * t81) * t53 * t79;
t46 = (t66 * t49 - (-t69 * t82 + t55 * t65 + (-t55 * t77 + t82) * t52) * t65 * t50) * t71 * t48;
t1 = [t63 * t61 * t79, 0, t52, t52, 0, 0; (-t49 * t77 - (-t55 * t62 * t63 * t78 + (t61 - 0.1e1) * t65 * t54) * t65 * t83) * t48, 0, t46, t46, 0, 0; ((-t66 * t76 - t73) * t56 - (-t66 * t75 + t74) * t81) * t53, 0, t47, t47, t72 * t53, 0;];
Ja_rot  = t1;
