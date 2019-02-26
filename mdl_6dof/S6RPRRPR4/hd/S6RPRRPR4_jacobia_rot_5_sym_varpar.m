% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:47
% EndTime: 2019-02-26 21:02:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (530->21), mult. (305->53), div. (74->9), fcn. (454->9), ass. (0->37)
t63 = pkin(10) + qJ(3) + qJ(4);
t62 = cos(t63);
t61 = sin(t63);
t66 = sin(qJ(1));
t72 = t66 * t61;
t56 = atan2(-t72, -t62);
t50 = sin(t56);
t51 = cos(t56);
t47 = -t50 * t72 - t51 * t62;
t46 = 0.1e1 / t47 ^ 2;
t67 = cos(qJ(1));
t78 = t46 * t67 ^ 2;
t77 = t50 * t62;
t65 = cos(pkin(11));
t68 = t67 * t65;
t64 = sin(pkin(11));
t71 = t66 * t64;
t55 = t62 * t68 + t71;
t53 = 0.1e1 / t55 ^ 2;
t69 = t67 * t64;
t70 = t66 * t65;
t54 = t62 * t69 - t70;
t76 = t53 * t54;
t58 = t61 ^ 2;
t75 = t58 / t62 ^ 2;
t74 = t61 * t67;
t57 = 0.1e1 / (t66 ^ 2 * t75 + 0.1e1);
t73 = t66 * t57;
t59 = 0.1e1 / t62;
t52 = 0.1e1 / t55;
t49 = 0.1e1 / (t54 ^ 2 * t53 + 0.1e1);
t48 = (0.1e1 + t75) * t73;
t45 = 0.1e1 / t47;
t44 = 0.1e1 / (t58 * t78 + 0.1e1);
t43 = (-t52 * t64 + t65 * t76) * t49 * t74;
t42 = (t62 * t45 - (-t66 * t77 + t51 * t61 + (-t51 * t72 + t77) * t48) * t61 * t46) * t67 * t44;
t1 = [t59 * t57 * t74, 0, t48, t48, 0, 0; (-t45 * t72 - (-t51 * t58 * t59 * t73 + (t57 - 0.1e1) * t61 * t50) * t61 * t78) * t44, 0, t42, t42, 0, 0; ((-t62 * t71 - t68) * t52 - (-t62 * t70 + t69) * t76) * t49, 0, t43, t43, 0, 0;];
Ja_rot  = t1;
