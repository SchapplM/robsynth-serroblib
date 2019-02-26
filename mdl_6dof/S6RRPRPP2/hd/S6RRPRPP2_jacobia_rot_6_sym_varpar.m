% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPP2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:35:21
% EndTime: 2019-02-26 21:35:21
% DurationCPUTime: 0.10s
% Computational Cost: add. (194->20), mult. (227->54), div. (53->9), fcn. (337->9), ass. (0->36)
t51 = qJ(2) + pkin(9);
t50 = cos(t51);
t49 = sin(t51);
t53 = sin(qJ(1));
t60 = t53 * t49;
t45 = atan2(t60, t50);
t42 = sin(t45);
t43 = cos(t45);
t34 = t42 * t60 + t43 * t50;
t33 = 0.1e1 / t34 ^ 2;
t55 = cos(qJ(1));
t67 = t33 * t55 ^ 2;
t54 = cos(qJ(4));
t56 = t55 * t54;
t52 = sin(qJ(4));
t59 = t53 * t52;
t41 = t50 * t56 + t59;
t39 = 0.1e1 / t41 ^ 2;
t57 = t55 * t52;
t58 = t53 * t54;
t40 = -t50 * t57 + t58;
t66 = t40 ^ 2 * t39;
t65 = t39 * t40;
t64 = t42 * t50;
t46 = t49 ^ 2;
t63 = t46 / t50 ^ 2;
t62 = t49 * t55;
t44 = 0.1e1 / (t53 ^ 2 * t63 + 0.1e1);
t61 = t53 * t44;
t47 = 0.1e1 / t50;
t38 = 0.1e1 / t41;
t36 = 0.1e1 / (0.1e1 + t66);
t35 = (0.1e1 + t63) * t61;
t32 = 0.1e1 / t34;
t31 = 0.1e1 / (t46 * t67 + 0.1e1);
t1 = [t47 * t44 * t62, t35, 0, 0, 0, 0; (t32 * t60 + (t43 * t46 * t47 * t61 + (-t44 + 0.1e1) * t49 * t42) * t49 * t67) * t31 (-t50 * t32 + (t53 * t64 - t43 * t49 + (t43 * t60 - t64) * t35) * t49 * t33) * t55 * t31, 0, 0, 0, 0; ((t50 * t59 + t56) * t38 - (-t50 * t58 + t57) * t65) * t36 (t38 * t52 + t54 * t65) * t36 * t62, 0 (-t38 * t41 - t66) * t36, 0, 0;];
Ja_rot  = t1;
