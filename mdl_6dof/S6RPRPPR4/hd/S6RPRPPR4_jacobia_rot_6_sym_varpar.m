% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:53
% EndTime: 2019-02-26 20:40:54
% DurationCPUTime: 0.12s
% Computational Cost: add. (265->25), mult. (347->66), div. (52->9), fcn. (503->11), ass. (0->42)
t55 = pkin(9) + qJ(3);
t54 = cos(t55);
t53 = sin(t55);
t59 = sin(qJ(1));
t67 = t59 * t53;
t49 = atan2(t67, t54);
t46 = sin(t49);
t47 = cos(t49);
t36 = t46 * t67 + t47 * t54;
t35 = 0.1e1 / t36 ^ 2;
t61 = cos(qJ(1));
t73 = t35 * t61 ^ 2;
t56 = sin(pkin(10));
t64 = t61 * t56;
t57 = cos(pkin(10));
t65 = t59 * t57;
t44 = t54 * t64 - t65;
t63 = t61 * t57;
t66 = t59 * t56;
t45 = t54 * t63 + t66;
t58 = sin(qJ(6));
t60 = cos(qJ(6));
t41 = t44 * t58 + t45 * t60;
t39 = 0.1e1 / t41 ^ 2;
t40 = -t44 * t60 + t45 * t58;
t72 = t39 * t40;
t71 = t46 * t54;
t50 = t53 ^ 2;
t70 = t50 / t54 ^ 2;
t69 = t53 * t61;
t48 = 0.1e1 / (t59 ^ 2 * t70 + 0.1e1);
t68 = t59 * t48;
t62 = t40 ^ 2 * t39 + 0.1e1;
t51 = 0.1e1 / t54;
t43 = -t54 * t65 + t64;
t42 = -t54 * t66 - t63;
t38 = 0.1e1 / t41;
t37 = (0.1e1 + t70) * t68;
t34 = 0.1e1 / t36;
t33 = 0.1e1 / (t50 * t73 + 0.1e1);
t32 = 0.1e1 / t62;
t1 = [t51 * t48 * t69, 0, t37, 0, 0, 0; (t34 * t67 + (t47 * t50 * t51 * t68 + (-t48 + 0.1e1) * t53 * t46) * t53 * t73) * t33, 0 (-t54 * t34 + (t59 * t71 - t47 * t53 + (t47 * t67 - t71) * t37) * t53 * t35) * t61 * t33, 0, 0, 0; ((-t42 * t60 + t43 * t58) * t38 - (t42 * t58 + t43 * t60) * t72) * t32, 0 ((t56 * t60 - t57 * t58) * t38 - (-t56 * t58 - t57 * t60) * t72) * t32 * t69, 0, 0, t62 * t32;];
Ja_rot  = t1;
