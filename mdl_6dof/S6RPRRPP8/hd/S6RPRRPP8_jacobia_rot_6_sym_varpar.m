% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_jacobia_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:00:20
% EndTime: 2019-02-26 21:00:20
% DurationCPUTime: 0.09s
% Computational Cost: add. (157->24), mult. (486->69), div. (108->11), fcn. (755->9), ass. (0->39)
t50 = sin(qJ(3));
t52 = cos(qJ(4));
t54 = cos(qJ(1));
t55 = t54 * t52;
t49 = sin(qJ(4));
t51 = sin(qJ(1));
t59 = t51 * t49;
t41 = t50 * t55 - t59;
t53 = cos(qJ(3));
t57 = t53 * t52;
t37 = atan2(t41, t57);
t33 = sin(t37);
t34 = cos(t37);
t32 = t33 * t41 + t34 * t57;
t31 = 0.1e1 / t32 ^ 2;
t56 = t54 * t49;
t58 = t51 * t52;
t39 = t50 * t58 + t56;
t66 = t31 * t39;
t64 = t34 * t41;
t38 = -t50 * t59 + t55;
t44 = 0.1e1 / t51 ^ 2;
t48 = 0.1e1 / t53 ^ 2;
t35 = 0.1e1 / (t38 ^ 2 * t44 * t48 + 0.1e1);
t47 = 0.1e1 / t53;
t63 = t35 * t47;
t62 = t39 ^ 2 * t31;
t45 = 0.1e1 / t52;
t61 = t45 * t47;
t60 = t48 * t50;
t46 = 0.1e1 / t52 ^ 2;
t43 = 0.1e1 / t51;
t40 = t50 * t56 + t58;
t36 = 0.1e1 / (t41 ^ 2 * t48 * t46 + 0.1e1);
t30 = 0.1e1 / t32;
t29 = (t41 * t45 * t60 + t54) * t36;
t28 = 0.1e1 / (0.1e1 + t62);
t27 = (t41 * t46 * t49 - t40 * t45) * t47 * t36;
t1 = [-t39 * t36 * t61, 0, t29, t27, 0, 0; (t41 * t30 - (-t33 + (-t61 * t64 + t33) * t36) * t62) * t28, 0 (-t29 * t64 * t66 + (t51 * t53 * t30 - (-t34 * t50 + (-t29 + t54) * t53 * t33) * t66) * t52) * t28 (t38 * t30 - (-t34 * t53 * t49 - t33 * t40 + (-t33 * t57 + t64) * t27) * t66) * t28, 0, 0; (t38 * t44 * t54 + t40 * t43) * t63, 0 (-t38 * t43 * t60 + t49) * t35, t39 * t43 * t63, 0, 0;];
Ja_rot  = t1;
