% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR7
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR7_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobia_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:18
% EndTime: 2019-02-26 22:19:18
% DurationCPUTime: 0.11s
% Computational Cost: add. (221->25), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->41)
t59 = sin(qJ(2));
t61 = cos(qJ(2));
t62 = cos(qJ(1));
t60 = sin(qJ(1));
t66 = cos(pkin(6));
t64 = t60 * t66;
t50 = -t59 * t64 + t62 * t61;
t55 = qJ(3) + pkin(12);
t52 = sin(t55);
t53 = cos(t55);
t58 = sin(pkin(6));
t69 = t58 * t60;
t41 = t50 * t53 + t52 * t69;
t39 = 0.1e1 / t41 ^ 2;
t40 = t50 * t52 - t53 * t69;
t73 = t39 * t40;
t63 = t62 * t66;
t46 = t60 * t59 - t61 * t63;
t68 = t58 * t61;
t44 = atan2(-t46, -t68);
t43 = cos(t44);
t72 = t43 * t46;
t42 = sin(t44);
t36 = -t42 * t46 - t43 * t68;
t35 = 0.1e1 / t36 ^ 2;
t49 = t62 * t59 + t61 * t64;
t71 = t49 ^ 2 * t35;
t54 = 0.1e1 / t58;
t56 = 0.1e1 / t61;
t70 = t54 * t56;
t67 = t58 * t62;
t65 = t40 ^ 2 * t39 + 0.1e1;
t57 = 0.1e1 / t61 ^ 2;
t48 = t59 * t63 + t60 * t61;
t45 = 0.1e1 / (0.1e1 + t46 ^ 2 / t58 ^ 2 * t57);
t38 = 0.1e1 / t41;
t37 = 0.1e1 / t65;
t34 = 0.1e1 / t36;
t33 = 0.1e1 / (0.1e1 + t71);
t32 = (t46 * t57 * t59 + t48 * t56) * t54 * t45;
t1 = [t49 * t45 * t70, t32, 0, 0, 0, 0; (-t46 * t34 - (-t42 + (-t70 * t72 + t42) * t45) * t71) * t33 (t50 * t34 - (t43 * t58 * t59 - t42 * t48 + (t42 * t68 - t72) * t32) * t49 * t35) * t33, 0, 0, 0, 0; ((-t48 * t52 - t53 * t67) * t38 - (-t48 * t53 + t52 * t67) * t73) * t37 (-t52 * t38 + t53 * t73) * t49 * t37, t65 * t37, 0, 0, 0;];
Ja_rot  = t1;
