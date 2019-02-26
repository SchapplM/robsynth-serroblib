% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP3
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

function Ja_rot = S6RPRRRP3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:09:12
% EndTime: 2019-02-26 21:09:12
% DurationCPUTime: 0.10s
% Computational Cost: add. (304->22), mult. (251->57), div. (57->9), fcn. (367->9), ass. (0->37)
t54 = qJ(1) + pkin(10);
t51 = cos(t54);
t71 = t51 ^ 2;
t60 = cos(qJ(3));
t50 = sin(t54);
t59 = sin(qJ(3));
t66 = t50 * t59;
t48 = atan2(-t66, -t60);
t46 = sin(t48);
t47 = cos(t48);
t40 = -t46 * t66 - t47 * t60;
t39 = 0.1e1 / t40 ^ 2;
t70 = t39 * t59;
t58 = qJ(4) + qJ(5);
t52 = sin(t58);
t53 = cos(t58);
t63 = t53 * t60;
t45 = t50 * t52 + t51 * t63;
t43 = 0.1e1 / t45 ^ 2;
t64 = t52 * t60;
t44 = -t50 * t53 + t51 * t64;
t69 = t43 * t44;
t68 = t46 * t60;
t55 = t59 ^ 2;
t62 = t55 / t60 ^ 2;
t49 = 0.1e1 / (t50 ^ 2 * t62 + 0.1e1);
t67 = t50 * t49;
t65 = t51 * t59;
t61 = t44 ^ 2 * t43 + 0.1e1;
t56 = 0.1e1 / t60;
t42 = 0.1e1 / t45;
t41 = (0.1e1 + t62) * t67;
t38 = 0.1e1 / t40;
t37 = 0.1e1 / t61;
t36 = 0.1e1 / (t71 * t55 * t39 + 0.1e1);
t35 = t61 * t37;
t1 = [t56 * t49 * t65, 0, t41, 0, 0, 0; (-t38 * t66 - (-t47 * t55 * t56 * t67 + (t49 - 0.1e1) * t59 * t46) * t71 * t70) * t36, 0 (t60 * t38 - (-t50 * t68 + t47 * t59 + (-t47 * t66 + t68) * t41) * t70) * t51 * t36, 0, 0, 0; ((-t50 * t64 - t51 * t53) * t42 - (-t50 * t63 + t51 * t52) * t69) * t37, 0 (-t42 * t52 + t53 * t69) * t37 * t65, t35, t35, 0;];
Ja_rot  = t1;
