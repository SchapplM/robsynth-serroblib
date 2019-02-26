% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_rot = S6RPRPPR4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:54
% EndTime: 2019-02-26 20:40:54
% DurationCPUTime: 0.08s
% Computational Cost: add. (227->21), mult. (329->51), div. (61->11), fcn. (494->9), ass. (0->36)
t49 = pkin(9) + qJ(3);
t46 = sin(t49);
t66 = t46 ^ 2;
t47 = cos(t49);
t52 = cos(pkin(10));
t54 = cos(qJ(1));
t56 = t54 * t52;
t51 = sin(pkin(10));
t53 = sin(qJ(1));
t59 = t53 * t51;
t38 = t47 * t59 + t56;
t60 = t46 * t51;
t34 = atan2(-t38, t60);
t31 = sin(t34);
t32 = cos(t34);
t30 = -t31 * t38 + t32 * t60;
t29 = 0.1e1 / t30 ^ 2;
t57 = t54 * t51;
t58 = t53 * t52;
t40 = t47 * t57 - t58;
t65 = t29 * t40;
t63 = t32 * t38;
t62 = t40 ^ 2 * t29;
t48 = 0.1e1 / t51;
t61 = 0.1e1 / t46 * t48;
t41 = t47 * t56 + t59;
t37 = 0.1e1 / t41 ^ 2;
t55 = t54 ^ 2 * t66 * t37;
t45 = 0.1e1 / t66;
t36 = 0.1e1 / t41;
t35 = 0.1e1 / (0.1e1 + t55);
t33 = 0.1e1 / (0.1e1 + t38 ^ 2 * t45 / t51 ^ 2);
t28 = 0.1e1 / t30;
t27 = (t38 * t45 * t47 * t48 + t53) * t33;
t26 = 0.1e1 / (0.1e1 + t62);
t1 = [-t40 * t33 * t61, 0, t27, 0, 0, 0; (-t38 * t28 - (-t31 + (t61 * t63 + t31) * t33) * t62) * t26, 0 (t27 * t63 * t65 + (-t54 * t46 * t28 - (t32 * t47 + (-t27 + t53) * t46 * t31) * t65) * t51) * t26, 0, 0, 0; (t53 * t36 + (-t47 * t58 + t57) * t54 * t37) * t46 * t35, 0 (-t36 * t47 * t54 - t52 * t55) * t35, 0, 0, 0;];
Ja_rot  = t1;
