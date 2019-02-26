% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP7
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
% Datum: 2019-02-26 21:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP7_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:11:19
% EndTime: 2019-02-26 21:11:19
% DurationCPUTime: 0.11s
% Computational Cost: add. (317->22), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->38)
t61 = pkin(10) + qJ(3);
t58 = cos(t61);
t57 = sin(t61);
t63 = sin(qJ(1));
t70 = t63 * t57;
t52 = atan2(-t70, -t58);
t50 = sin(t52);
t51 = cos(t52);
t43 = -t50 * t70 - t51 * t58;
t42 = 0.1e1 / t43 ^ 2;
t64 = cos(qJ(1));
t76 = t42 * t64 ^ 2;
t62 = qJ(4) + qJ(5);
t60 = cos(t62);
t66 = t64 * t60;
t59 = sin(t62);
t69 = t63 * t59;
t49 = t58 * t66 + t69;
t47 = 0.1e1 / t49 ^ 2;
t67 = t64 * t59;
t68 = t63 * t60;
t48 = t58 * t67 - t68;
t75 = t47 * t48;
t74 = t50 * t58;
t54 = t57 ^ 2;
t73 = t54 / t58 ^ 2;
t72 = t57 * t64;
t53 = 0.1e1 / (t63 ^ 2 * t73 + 0.1e1);
t71 = t63 * t53;
t65 = t48 ^ 2 * t47 + 0.1e1;
t55 = 0.1e1 / t58;
t46 = 0.1e1 / t49;
t45 = 0.1e1 / t65;
t44 = (0.1e1 + t73) * t71;
t41 = 0.1e1 / t43;
t40 = 0.1e1 / (t54 * t76 + 0.1e1);
t39 = t65 * t45;
t1 = [t55 * t53 * t72, 0, t44, 0, 0, 0; (-t41 * t70 - (-t51 * t54 * t55 * t71 + (t53 - 0.1e1) * t57 * t50) * t57 * t76) * t40, 0 (t58 * t41 - (-t63 * t74 + t51 * t57 + (-t51 * t70 + t74) * t44) * t57 * t42) * t64 * t40, 0, 0, 0; ((-t58 * t69 - t66) * t46 - (-t58 * t68 + t67) * t75) * t45, 0 (-t46 * t59 + t60 * t75) * t45 * t72, t39, t39, 0;];
Ja_rot  = t1;
