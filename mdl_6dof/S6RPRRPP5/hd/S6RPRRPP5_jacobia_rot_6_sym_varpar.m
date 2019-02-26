% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:58:38
% EndTime: 2019-02-26 20:58:38
% DurationCPUTime: 0.07s
% Computational Cost: add. (194->20), mult. (227->54), div. (53->9), fcn. (337->9), ass. (0->36)
t48 = pkin(9) + qJ(3);
t47 = cos(t48);
t46 = sin(t48);
t50 = sin(qJ(1));
t57 = t50 * t46;
t42 = atan2(t57, t47);
t39 = sin(t42);
t40 = cos(t42);
t31 = t39 * t57 + t40 * t47;
t30 = 0.1e1 / t31 ^ 2;
t52 = cos(qJ(1));
t64 = t30 * t52 ^ 2;
t51 = cos(qJ(4));
t53 = t52 * t51;
t49 = sin(qJ(4));
t56 = t50 * t49;
t38 = t47 * t53 + t56;
t36 = 0.1e1 / t38 ^ 2;
t54 = t52 * t49;
t55 = t50 * t51;
t37 = -t47 * t54 + t55;
t63 = t37 ^ 2 * t36;
t62 = t36 * t37;
t61 = t39 * t47;
t43 = t46 ^ 2;
t60 = t43 / t47 ^ 2;
t59 = t46 * t52;
t41 = 0.1e1 / (t50 ^ 2 * t60 + 0.1e1);
t58 = t50 * t41;
t44 = 0.1e1 / t47;
t35 = 0.1e1 / t38;
t33 = 0.1e1 / (0.1e1 + t63);
t32 = (0.1e1 + t60) * t58;
t29 = 0.1e1 / t31;
t28 = 0.1e1 / (t43 * t64 + 0.1e1);
t1 = [t44 * t41 * t59, 0, t32, 0, 0, 0; (t29 * t57 + (t40 * t43 * t44 * t58 + (-t41 + 0.1e1) * t46 * t39) * t46 * t64) * t28, 0 (-t47 * t29 + (t50 * t61 - t40 * t46 + (t40 * t57 - t61) * t32) * t46 * t30) * t52 * t28, 0, 0, 0; ((t47 * t56 + t53) * t35 - (-t47 * t55 + t54) * t62) * t33, 0 (t35 * t49 + t51 * t62) * t33 * t59 (-t35 * t38 - t63) * t33, 0, 0;];
Ja_rot  = t1;
