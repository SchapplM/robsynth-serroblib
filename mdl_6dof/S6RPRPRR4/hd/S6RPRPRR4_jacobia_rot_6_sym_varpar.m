% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:50:49
% EndTime: 2019-02-26 20:50:49
% DurationCPUTime: 0.10s
% Computational Cost: add. (282->22), mult. (251->57), div. (57->9), fcn. (367->9), ass. (0->37)
t56 = qJ(1) + pkin(10);
t53 = cos(t56);
t73 = t53 ^ 2;
t61 = sin(qJ(3));
t52 = sin(t56);
t62 = cos(qJ(3));
t68 = t52 * t62;
t50 = atan2(-t68, t61);
t48 = sin(t50);
t49 = cos(t50);
t43 = -t48 * t68 + t49 * t61;
t41 = 0.1e1 / t43 ^ 2;
t72 = t41 * t62;
t60 = qJ(5) + qJ(6);
t55 = cos(t60);
t54 = sin(t60);
t66 = t54 * t61;
t47 = t52 * t55 + t53 * t66;
t45 = 0.1e1 / t47 ^ 2;
t65 = t55 * t61;
t46 = t52 * t54 - t53 * t65;
t71 = t45 * t46;
t70 = t48 * t61;
t59 = t62 ^ 2;
t64 = 0.1e1 / t61 ^ 2 * t59;
t51 = 0.1e1 / (t52 ^ 2 * t64 + 0.1e1);
t69 = t52 * t51;
t67 = t53 * t62;
t63 = t46 ^ 2 * t45 + 0.1e1;
t57 = 0.1e1 / t61;
t44 = 0.1e1 / t47;
t42 = (0.1e1 + t64) * t69;
t40 = 0.1e1 / t43;
t39 = 0.1e1 / t63;
t38 = 0.1e1 / (t73 * t59 * t41 + 0.1e1);
t37 = t63 * t39;
t1 = [-t57 * t51 * t67, 0, t42, 0, 0, 0; (-t40 * t68 - (t49 * t57 * t59 * t69 + (t51 - 0.1e1) * t62 * t48) * t73 * t72) * t38, 0 (-t61 * t40 - (t52 * t70 + t49 * t62 + (-t49 * t68 - t70) * t42) * t72) * t53 * t38, 0, 0, 0; ((t52 * t65 + t53 * t54) * t44 - (-t52 * t66 + t53 * t55) * t71) * t39, 0 (-t44 * t55 - t54 * t71) * t39 * t67, 0, t37, t37;];
Ja_rot  = t1;
