% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR3
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
% Datum: 2019-02-26 20:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:50:14
% EndTime: 2019-02-26 20:50:14
% DurationCPUTime: 0.11s
% Computational Cost: add. (366->22), mult. (251->57), div. (57->9), fcn. (367->9), ass. (0->36)
t56 = qJ(1) + pkin(10);
t54 = cos(t56);
t71 = t54 ^ 2;
t61 = cos(qJ(3));
t53 = sin(t56);
t60 = sin(qJ(3));
t67 = t53 * t60;
t49 = atan2(-t67, -t61);
t47 = sin(t49);
t48 = cos(t49);
t41 = -t47 * t67 - t48 * t61;
t40 = 0.1e1 / t41 ^ 2;
t70 = t40 * t60;
t55 = pkin(11) + qJ(5) + qJ(6);
t51 = sin(t55);
t52 = cos(t55);
t64 = t54 * t61;
t46 = t53 * t51 + t52 * t64;
t44 = 0.1e1 / t46 ^ 2;
t45 = t51 * t64 - t53 * t52;
t69 = t44 * t45;
t57 = t60 ^ 2;
t63 = t57 / t61 ^ 2;
t50 = 0.1e1 / (t53 ^ 2 * t63 + 0.1e1);
t68 = t53 * t50;
t66 = t53 * t61;
t65 = t54 * t60;
t62 = t45 ^ 2 * t44 + 0.1e1;
t58 = 0.1e1 / t61;
t43 = 0.1e1 / t46;
t42 = (0.1e1 + t63) * t68;
t39 = 0.1e1 / t41;
t38 = 0.1e1 / (t71 * t57 * t40 + 0.1e1);
t37 = 0.1e1 / t62;
t36 = t62 * t37;
t1 = [t58 * t50 * t65, 0, t42, 0, 0, 0; (-t39 * t67 - (-t48 * t57 * t58 * t68 + (t50 - 0.1e1) * t60 * t47) * t71 * t70) * t38, 0 (t61 * t39 - (-t47 * t66 + t48 * t60 + (t47 * t61 - t48 * t67) * t42) * t70) * t54 * t38, 0, 0, 0; ((-t51 * t66 - t54 * t52) * t43 - (t54 * t51 - t52 * t66) * t69) * t37, 0 (-t43 * t51 + t52 * t69) * t37 * t65, 0, t36, t36;];
Ja_rot  = t1;
