% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:51:42
% EndTime: 2019-02-26 20:51:42
% DurationCPUTime: 0.08s
% Computational Cost: add. (379->22), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->38)
t64 = pkin(10) + qJ(3);
t62 = cos(t64);
t61 = sin(t64);
t65 = sin(qJ(1));
t70 = t65 * t61;
t54 = atan2(-t70, -t62);
t52 = sin(t54);
t53 = cos(t54);
t46 = -t52 * t70 - t53 * t62;
t45 = 0.1e1 / t46 ^ 2;
t66 = cos(qJ(1));
t78 = t45 * t66 ^ 2;
t63 = pkin(11) + qJ(5) + qJ(6);
t57 = cos(t63);
t68 = t66 * t57;
t56 = sin(t63);
t72 = t65 * t56;
t51 = t62 * t68 + t72;
t49 = 0.1e1 / t51 ^ 2;
t69 = t66 * t56;
t71 = t65 * t57;
t50 = t62 * t69 - t71;
t77 = t49 * t50;
t76 = t52 * t62;
t58 = t61 ^ 2;
t75 = t58 / t62 ^ 2;
t74 = t61 * t66;
t55 = 0.1e1 / (t65 ^ 2 * t75 + 0.1e1);
t73 = t65 * t55;
t67 = t50 ^ 2 * t49 + 0.1e1;
t59 = 0.1e1 / t62;
t48 = 0.1e1 / t51;
t47 = (0.1e1 + t75) * t73;
t44 = 0.1e1 / t46;
t43 = 0.1e1 / t67;
t42 = 0.1e1 / (t58 * t78 + 0.1e1);
t41 = t67 * t43;
t1 = [t59 * t55 * t74, 0, t47, 0, 0, 0; (-t44 * t70 - (-t53 * t58 * t59 * t73 + (t55 - 0.1e1) * t61 * t52) * t61 * t78) * t42, 0 (t62 * t44 - (-t65 * t76 + t53 * t61 + (-t53 * t70 + t76) * t47) * t61 * t45) * t66 * t42, 0, 0, 0; ((-t62 * t72 - t68) * t48 - (-t62 * t71 + t69) * t77) * t43, 0 (-t48 * t56 + t57 * t77) * t43 * t74, 0, t41, t41;];
Ja_rot  = t1;
