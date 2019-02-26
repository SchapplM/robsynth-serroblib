% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRP8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:34:12
% EndTime: 2019-02-26 20:34:12
% DurationCPUTime: 0.10s
% Computational Cost: add. (354->25), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->41)
t52 = pkin(9) + qJ(4);
t51 = cos(t52);
t72 = t51 ^ 2;
t50 = sin(t52);
t56 = sin(qJ(5));
t59 = cos(qJ(1));
t62 = t59 * t56;
t57 = sin(qJ(1));
t58 = cos(qJ(5));
t63 = t57 * t58;
t44 = t50 * t62 + t63;
t66 = t51 * t56;
t39 = atan2(t44, t66);
t35 = sin(t39);
t36 = cos(t39);
t34 = t35 * t44 + t36 * t66;
t33 = 0.1e1 / t34 ^ 2;
t61 = t59 * t58;
t64 = t57 * t56;
t42 = t50 * t64 - t61;
t71 = t33 * t42;
t69 = t36 * t44;
t68 = t42 ^ 2 * t33;
t48 = 0.1e1 / t51;
t53 = 0.1e1 / t56;
t67 = t48 * t53;
t65 = t51 * t57;
t43 = t50 * t63 + t62;
t41 = 0.1e1 / t43 ^ 2;
t60 = t57 ^ 2 * t72 * t41;
t54 = 0.1e1 / t56 ^ 2;
t49 = 0.1e1 / t72;
t45 = t50 * t61 - t64;
t40 = 0.1e1 / t43;
t38 = 0.1e1 / (t44 ^ 2 * t49 * t54 + 0.1e1);
t37 = 0.1e1 / (0.1e1 + t60);
t32 = 0.1e1 / t34;
t31 = (t44 * t49 * t50 * t53 + t59) * t38;
t30 = 0.1e1 / (0.1e1 + t68);
t29 = (-t44 * t54 * t58 + t45 * t53) * t48 * t38;
t1 = [-t42 * t38 * t67, 0, 0, t31, t29, 0; (t44 * t32 - (-t35 + (-t67 * t69 + t35) * t38) * t68) * t30, 0, 0 (-t31 * t69 * t71 + (t32 * t65 - (-t36 * t50 + (-t31 + t59) * t51 * t35) * t71) * t56) * t30 (t43 * t32 - (t36 * t51 * t58 + t35 * t45 + (-t35 * t66 + t69) * t29) * t71) * t30, 0; (-t41 * t45 * t57 + t40 * t59) * t51 * t37, 0, 0 (-t40 * t50 * t57 - t58 * t60) * t37, t42 * t41 * t37 * t65, 0;];
Ja_rot  = t1;
