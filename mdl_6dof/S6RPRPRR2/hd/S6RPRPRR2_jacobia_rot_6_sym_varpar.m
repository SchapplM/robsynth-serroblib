% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR2
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
% Datum: 2019-02-26 20:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:49:38
% EndTime: 2019-02-26 20:49:39
% DurationCPUTime: 0.11s
% Computational Cost: add. (440->23), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->39)
t67 = qJ(3) + pkin(11);
t63 = cos(t67);
t61 = sin(t67);
t68 = qJ(1) + pkin(10);
t62 = sin(t68);
t75 = t62 * t61;
t56 = atan2(-t75, -t63);
t54 = sin(t56);
t55 = cos(t56);
t47 = -t54 * t75 - t55 * t63;
t46 = 0.1e1 / t47 ^ 2;
t64 = cos(t68);
t81 = t46 * t64 ^ 2;
t69 = qJ(5) + qJ(6);
t66 = cos(t69);
t71 = t64 * t66;
t65 = sin(t69);
t74 = t62 * t65;
t53 = t63 * t71 + t74;
t51 = 0.1e1 / t53 ^ 2;
t72 = t64 * t65;
t73 = t62 * t66;
t52 = t63 * t72 - t73;
t80 = t51 * t52;
t79 = t54 * t63;
t58 = t61 ^ 2;
t78 = t58 / t63 ^ 2;
t77 = t61 * t64;
t57 = 0.1e1 / (t62 ^ 2 * t78 + 0.1e1);
t76 = t62 * t57;
t70 = t52 ^ 2 * t51 + 0.1e1;
t59 = 0.1e1 / t63;
t50 = 0.1e1 / t53;
t49 = 0.1e1 / t70;
t48 = (0.1e1 + t78) * t76;
t45 = 0.1e1 / t47;
t44 = 0.1e1 / (t58 * t81 + 0.1e1);
t43 = t70 * t49;
t1 = [t59 * t57 * t77, 0, t48, 0, 0, 0; (-t45 * t75 - (-t55 * t58 * t59 * t76 + (t57 - 0.1e1) * t61 * t54) * t61 * t81) * t44, 0 (t63 * t45 - (-t62 * t79 + t55 * t61 + (-t55 * t75 + t79) * t48) * t61 * t46) * t64 * t44, 0, 0, 0; ((-t63 * t74 - t71) * t50 - (-t63 * t73 + t72) * t80) * t49, 0 (-t50 * t65 + t66 * t80) * t49 * t77, 0, t43, t43;];
Ja_rot  = t1;
