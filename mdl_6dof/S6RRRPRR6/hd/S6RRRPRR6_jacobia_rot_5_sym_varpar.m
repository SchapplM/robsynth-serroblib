% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR6_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:47
% EndTime: 2019-02-26 22:18:48
% DurationCPUTime: 0.11s
% Computational Cost: add. (243->21), mult. (251->55), div. (57->9), fcn. (367->9), ass. (0->35)
t58 = cos(qJ(2));
t56 = sin(qJ(2));
t57 = sin(qJ(1));
t64 = t57 * t56;
t48 = atan2(-t64, -t58);
t46 = sin(t48);
t47 = cos(t48);
t40 = -t46 * t64 - t47 * t58;
t39 = 0.1e1 / t40 ^ 2;
t59 = cos(qJ(1));
t69 = t39 * t59 ^ 2;
t52 = qJ(3) + pkin(11) + qJ(5);
t50 = sin(t52);
t51 = cos(t52);
t61 = t59 * t51;
t45 = t57 * t50 + t58 * t61;
t43 = 0.1e1 / t45 ^ 2;
t62 = t59 * t50;
t44 = -t57 * t51 + t58 * t62;
t68 = t43 * t44;
t53 = t56 ^ 2;
t67 = t53 / t58 ^ 2;
t66 = t56 * t59;
t49 = 0.1e1 / (t57 ^ 2 * t67 + 0.1e1);
t65 = t57 * t49;
t63 = t57 * t58;
t60 = t44 ^ 2 * t43 + 0.1e1;
t54 = 0.1e1 / t58;
t42 = 0.1e1 / t45;
t41 = (0.1e1 + t67) * t65;
t38 = 0.1e1 / t40;
t37 = 0.1e1 / (t53 * t69 + 0.1e1);
t36 = 0.1e1 / t60;
t35 = t60 * t36;
t1 = [t54 * t49 * t66, t41, 0, 0, 0, 0; (-t38 * t64 - (-t47 * t53 * t54 * t65 + (t49 - 0.1e1) * t56 * t46) * t56 * t69) * t37 (t58 * t38 - (-t46 * t63 + t47 * t56 + (t46 * t58 - t47 * t64) * t41) * t56 * t39) * t59 * t37, 0, 0, 0, 0; ((-t50 * t63 - t61) * t42 - (-t51 * t63 + t62) * t68) * t36 (-t42 * t50 + t51 * t68) * t36 * t66, t35, 0, t35, 0;];
Ja_rot  = t1;
