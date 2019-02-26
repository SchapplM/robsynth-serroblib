% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPPR3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:59:23
% EndTime: 2019-02-26 19:59:23
% DurationCPUTime: 0.09s
% Computational Cost: add. (91->17), mult. (274->44), div. (48->11), fcn. (404->11), ass. (0->30)
t51 = sin(pkin(10));
t53 = cos(pkin(10));
t58 = cos(qJ(2));
t54 = cos(pkin(6));
t56 = sin(qJ(2));
t60 = t54 * t56;
t47 = -t51 * t60 + t53 * t58;
t55 = sin(qJ(3));
t57 = cos(qJ(3));
t52 = sin(pkin(6));
t62 = t51 * t52;
t38 = t47 * t55 - t57 * t62;
t36 = 0.1e1 / t38 ^ 2;
t39 = t47 * t57 + t55 * t62;
t63 = t39 ^ 2 * t36;
t61 = t52 * t58;
t59 = t54 * t58;
t50 = 0.1e1 / t58 ^ 2;
t46 = -t51 * t59 - t53 * t56;
t45 = t51 * t58 + t53 * t60;
t44 = t51 * t56 - t53 * t59;
t43 = atan2(t44, t61);
t41 = cos(t43);
t40 = sin(t43);
t35 = 0.1e1 / t38;
t34 = 0.1e1 / (0.1e1 + t63);
t33 = t40 * t44 + t41 * t61;
t32 = 0.1e1 / t33 ^ 2;
t30 = (t45 / t58 + t56 * t44 * t50) / t52 / (0.1e1 + t44 ^ 2 / t52 ^ 2 * t50);
t1 = [0, t30, 0, 0, 0, 0; 0 (-t47 / t33 - (-t41 * t52 * t56 + t40 * t45 + (-t40 * t61 + t41 * t44) * t30) * t46 * t32) / (t46 ^ 2 * t32 + 0.1e1) 0, 0, 0, 0; 0 (-t36 * t39 * t55 + t35 * t57) * t46 * t34 (-t35 * t38 - t63) * t34, 0, 0, 0;];
Ja_rot  = t1;
