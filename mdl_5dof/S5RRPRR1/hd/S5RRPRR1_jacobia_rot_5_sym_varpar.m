% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
%
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRPRR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_jacobia_rot_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_jacobia_rot_5_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:22:52
% EndTime: 2019-07-18 17:22:52
% DurationCPUTime: 0.12s
% Computational Cost: add. (358->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
t67 = qJ(2) + qJ(4);
t66 = cos(t67);
t65 = sin(t67);
t69 = sin(qJ(1));
t77 = t69 * t65;
t60 = atan2(-t77, -t66);
t58 = sin(t60);
t59 = cos(t60);
t51 = -t58 * t77 - t59 * t66;
t50 = 0.1e1 / t51 ^ 2;
t71 = cos(qJ(1));
t83 = t50 * t71 ^ 2;
t70 = cos(qJ(5));
t73 = t71 * t70;
t68 = sin(qJ(5));
t76 = t69 * t68;
t57 = t66 * t73 + t76;
t55 = 0.1e1 / t57 ^ 2;
t74 = t71 * t68;
t75 = t69 * t70;
t56 = t66 * t74 - t75;
t82 = t55 * t56;
t81 = t58 * t66;
t62 = t65 ^ 2;
t80 = t62 / t66 ^ 2;
t79 = t65 * t71;
t61 = 0.1e1 / (t69 ^ 2 * t80 + 0.1e1);
t78 = t69 * t61;
t72 = t56 ^ 2 * t55 + 0.1e1;
t63 = 0.1e1 / t66;
t54 = 0.1e1 / t57;
t53 = 0.1e1 / t72;
t52 = (0.1e1 + t80) * t78;
t49 = 0.1e1 / t51;
t48 = 0.1e1 / (t62 * t83 + 0.1e1);
t47 = (-t54 * t68 + t70 * t82) * t53 * t79;
t46 = (t66 * t49 - (-t69 * t81 + t59 * t65 + (-t59 * t77 + t81) * t52) * t65 * t50) * t71 * t48;
t1 = [t63 * t61 * t79, t52, 0, t52, 0; (-t49 * t77 - (-t59 * t62 * t63 * t78 + (t61 - 0.1e1) * t65 * t58) * t65 * t83) * t48, t46, 0, t46, 0; ((-t66 * t76 - t73) * t54 - (-t66 * t75 + t74) * t82) * t53, t47, 0, t47, t72 * t53;];
Ja_rot  = t1;
