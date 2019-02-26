% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_rot = S6RPRRRP4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:09:38
% EndTime: 2019-02-26 21:09:38
% DurationCPUTime: 0.11s
% Computational Cost: add. (554->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
t69 = pkin(10) + qJ(3) + qJ(4);
t68 = cos(t69);
t67 = sin(t69);
t71 = sin(qJ(1));
t79 = t71 * t67;
t60 = atan2(-t79, -t68);
t56 = sin(t60);
t57 = cos(t60);
t53 = -t56 * t79 - t57 * t68;
t52 = 0.1e1 / t53 ^ 2;
t73 = cos(qJ(1));
t85 = t52 * t73 ^ 2;
t84 = t56 * t68;
t72 = cos(qJ(5));
t75 = t73 * t72;
t70 = sin(qJ(5));
t78 = t71 * t70;
t62 = t68 * t75 + t78;
t59 = 0.1e1 / t62 ^ 2;
t76 = t73 * t70;
t77 = t71 * t72;
t61 = t68 * t76 - t77;
t83 = t59 * t61;
t64 = t67 ^ 2;
t82 = t64 / t68 ^ 2;
t81 = t67 * t73;
t63 = 0.1e1 / (t71 ^ 2 * t82 + 0.1e1);
t80 = t71 * t63;
t74 = t61 ^ 2 * t59 + 0.1e1;
t65 = 0.1e1 / t68;
t58 = 0.1e1 / t62;
t55 = 0.1e1 / t74;
t54 = (0.1e1 + t82) * t80;
t51 = 0.1e1 / t53;
t50 = 0.1e1 / (t64 * t85 + 0.1e1);
t49 = (-t58 * t70 + t72 * t83) * t55 * t81;
t48 = (t68 * t51 - (-t71 * t84 + t57 * t67 + (-t57 * t79 + t84) * t54) * t67 * t52) * t73 * t50;
t1 = [t65 * t63 * t81, 0, t54, t54, 0, 0; (-t51 * t79 - (-t57 * t64 * t65 * t80 + (t63 - 0.1e1) * t67 * t56) * t67 * t85) * t50, 0, t48, t48, 0, 0; ((-t68 * t78 - t75) * t58 - (-t68 * t77 + t76) * t83) * t55, 0, t49, t49, t74 * t55, 0;];
Ja_rot  = t1;
