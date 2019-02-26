% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:43:41
% EndTime: 2019-02-26 20:43:41
% DurationCPUTime: 0.11s
% Computational Cost: add. (347->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
t53 = qJ(3) + pkin(10);
t51 = cos(t53);
t49 = sin(t53);
t54 = qJ(1) + pkin(9);
t50 = sin(t54);
t62 = t50 * t49;
t44 = atan2(-t62, -t51);
t42 = sin(t44);
t43 = cos(t44);
t35 = -t42 * t62 - t43 * t51;
t34 = 0.1e1 / t35 ^ 2;
t52 = cos(t54);
t68 = t34 * t52 ^ 2;
t56 = cos(qJ(5));
t58 = t52 * t56;
t55 = sin(qJ(5));
t61 = t50 * t55;
t41 = t51 * t58 + t61;
t39 = 0.1e1 / t41 ^ 2;
t59 = t52 * t55;
t60 = t50 * t56;
t40 = t51 * t59 - t60;
t67 = t39 * t40;
t66 = t42 * t51;
t46 = t49 ^ 2;
t65 = t46 / t51 ^ 2;
t64 = t49 * t52;
t45 = 0.1e1 / (t50 ^ 2 * t65 + 0.1e1);
t63 = t50 * t45;
t57 = t40 ^ 2 * t39 + 0.1e1;
t47 = 0.1e1 / t51;
t38 = 0.1e1 / t41;
t37 = 0.1e1 / t57;
t36 = (0.1e1 + t65) * t63;
t33 = 0.1e1 / t35;
t32 = 0.1e1 / (t46 * t68 + 0.1e1);
t1 = [t47 * t45 * t64, 0, t36, 0, 0, 0; (-t33 * t62 - (-t43 * t46 * t47 * t63 + (t45 - 0.1e1) * t49 * t42) * t49 * t68) * t32, 0 (t51 * t33 - (-t50 * t66 + t43 * t49 + (-t43 * t62 + t66) * t36) * t49 * t34) * t52 * t32, 0, 0, 0; ((-t51 * t61 - t58) * t38 - (-t51 * t60 + t59) * t67) * t37, 0 (-t38 * t55 + t56 * t67) * t37 * t64, 0, t57 * t37, 0;];
Ja_rot  = t1;
