% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP6_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:59:25
% EndTime: 2019-02-26 20:59:25
% DurationCPUTime: 0.10s
% Computational Cost: add. (134->20), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->37)
t46 = sin(qJ(1));
t61 = t46 ^ 2;
t45 = sin(qJ(3));
t47 = cos(qJ(3));
t48 = cos(qJ(1));
t50 = t48 * t47;
t37 = atan2(-t50, t45);
t35 = sin(t37);
t36 = cos(t37);
t29 = -t35 * t50 + t36 * t45;
t28 = 0.1e1 / t29 ^ 2;
t60 = t28 * t47;
t41 = qJ(4) + pkin(9);
t39 = sin(t41);
t52 = t48 * t39;
t40 = cos(t41);
t55 = t46 * t40;
t34 = t45 * t55 + t52;
t32 = 0.1e1 / t34 ^ 2;
t51 = t48 * t40;
t56 = t46 * t39;
t33 = t45 * t56 - t51;
t59 = t32 * t33;
t58 = t35 * t45;
t44 = t47 ^ 2;
t57 = 0.1e1 / t45 ^ 2 * t44;
t54 = t46 * t47;
t38 = 0.1e1 / (t48 ^ 2 * t57 + 0.1e1);
t53 = t48 * t38;
t49 = t33 ^ 2 * t32 + 0.1e1;
t42 = 0.1e1 / t45;
t31 = 0.1e1 / t34;
t30 = (0.1e1 + t57) * t53;
t27 = 0.1e1 / t29;
t26 = 0.1e1 / t49;
t25 = 0.1e1 / (t61 * t44 * t28 + 0.1e1);
t1 = [t42 * t38 * t54, 0, t30, 0, 0, 0; (-t27 * t50 + (-t36 * t42 * t44 * t53 + (-t38 + 0.1e1) * t47 * t35) * t61 * t60) * t25, 0 (t45 * t27 + (t48 * t58 + t36 * t47 + (-t36 * t50 - t58) * t30) * t60) * t46 * t25, 0, 0, 0; ((t45 * t52 + t55) * t31 - (t45 * t51 - t56) * t59) * t26, 0 (t31 * t39 - t40 * t59) * t26 * t54, t49 * t26, 0, 0;];
Ja_rot  = t1;
