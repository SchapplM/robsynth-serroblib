% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP9_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:48:11
% EndTime: 2019-02-26 20:48:11
% DurationCPUTime: 0.10s
% Computational Cost: add. (447->25), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->41)
t57 = cos(qJ(3));
t71 = t57 ^ 2;
t55 = sin(qJ(3));
t50 = pkin(9) + qJ(5);
t48 = sin(t50);
t58 = cos(qJ(1));
t61 = t58 * t48;
t49 = cos(t50);
t56 = sin(qJ(1));
t64 = t56 * t49;
t43 = t55 * t61 + t64;
t62 = t57 * t48;
t37 = atan2(t43, t62);
t34 = sin(t37);
t35 = cos(t37);
t33 = t34 * t43 + t35 * t62;
t32 = 0.1e1 / t33 ^ 2;
t60 = t58 * t49;
t65 = t56 * t48;
t41 = t55 * t65 - t60;
t70 = t32 * t41;
t68 = t35 * t43;
t67 = t41 ^ 2 * t32;
t46 = 0.1e1 / t48;
t53 = 0.1e1 / t57;
t66 = t46 * t53;
t63 = t56 * t57;
t42 = t55 * t64 + t61;
t40 = 0.1e1 / t42 ^ 2;
t59 = t56 ^ 2 * t71 * t40;
t54 = 0.1e1 / t71;
t47 = 0.1e1 / t48 ^ 2;
t44 = t55 * t60 - t65;
t39 = 0.1e1 / t42;
t38 = 0.1e1 / (0.1e1 + t59);
t36 = 0.1e1 / (t43 ^ 2 * t54 * t47 + 0.1e1);
t31 = 0.1e1 / t33;
t30 = (t43 * t46 * t54 * t55 + t58) * t36;
t29 = 0.1e1 / (0.1e1 + t67);
t28 = (-t43 * t47 * t49 + t44 * t46) * t53 * t36;
t1 = [-t41 * t36 * t66, 0, t30, 0, t28, 0; (t43 * t31 - (-t34 + (-t66 * t68 + t34) * t36) * t67) * t29, 0 (-t30 * t68 * t70 + (t31 * t63 - (-t35 * t55 + (-t30 + t58) * t57 * t34) * t70) * t48) * t29, 0 (t42 * t31 - (t35 * t57 * t49 + t34 * t44 + (-t34 * t62 + t68) * t28) * t70) * t29, 0; (-t40 * t44 * t56 + t39 * t58) * t57 * t38, 0 (-t39 * t55 * t56 - t49 * t59) * t38, 0, t41 * t40 * t38 * t63, 0;];
Ja_rot  = t1;
