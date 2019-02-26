% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRP2_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:47
% EndTime: 2019-02-26 20:01:47
% DurationCPUTime: 0.06s
% Computational Cost: add. (131->19), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->30)
t50 = sin(pkin(10));
t51 = sin(pkin(6));
t60 = t50 * t51;
t55 = cos(qJ(2));
t59 = t51 * t55;
t53 = cos(pkin(6));
t54 = sin(qJ(2));
t58 = t53 * t54;
t57 = t53 * t55;
t52 = cos(pkin(10));
t43 = -t50 * t58 + t52 * t55;
t48 = qJ(3) + pkin(11);
t45 = sin(t48);
t46 = cos(t48);
t34 = t43 * t46 + t45 * t60;
t32 = 0.1e1 / t34 ^ 2;
t33 = t43 * t45 - t46 * t60;
t56 = t33 ^ 2 * t32 + 0.1e1;
t49 = 0.1e1 / t55 ^ 2;
t42 = t50 * t57 + t52 * t54;
t41 = t50 * t55 + t52 * t58;
t39 = t50 * t54 - t52 * t57;
t37 = atan2(-t39, -t59);
t36 = cos(t37);
t35 = sin(t37);
t31 = 0.1e1 / t56;
t30 = -t35 * t39 - t36 * t59;
t29 = 0.1e1 / t30 ^ 2;
t27 = (t41 / t55 + t54 * t39 * t49) / t51 / (0.1e1 + t39 ^ 2 / t51 ^ 2 * t49);
t1 = [0, t27, 0, 0, 0, 0; 0 (t43 / t30 - (t36 * t51 * t54 - t35 * t41 + (t35 * t59 - t36 * t39) * t27) * t42 * t29) / (t42 ^ 2 * t29 + 0.1e1) 0, 0, 0, 0; 0 (-t45 / t34 + t46 * t33 * t32) * t42 * t31, t56 * t31, 0, 0, 0;];
Ja_rot  = t1;
