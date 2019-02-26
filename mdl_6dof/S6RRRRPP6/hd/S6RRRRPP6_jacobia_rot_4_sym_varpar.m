% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP6_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_jacobia_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:28:14
% EndTime: 2019-02-26 22:28:14
% DurationCPUTime: 0.11s
% Computational Cost: add. (181->21), mult. (251->55), div. (57->9), fcn. (367->9), ass. (0->35)
t56 = cos(qJ(2));
t54 = sin(qJ(2));
t55 = sin(qJ(1));
t62 = t55 * t54;
t46 = atan2(-t62, -t56);
t44 = sin(t46);
t45 = cos(t46);
t38 = -t44 * t62 - t45 * t56;
t37 = 0.1e1 / t38 ^ 2;
t57 = cos(qJ(1));
t67 = t37 * t57 ^ 2;
t53 = qJ(3) + qJ(4);
t48 = sin(t53);
t49 = cos(t53);
t59 = t57 * t49;
t43 = t55 * t48 + t56 * t59;
t41 = 0.1e1 / t43 ^ 2;
t60 = t57 * t48;
t42 = -t55 * t49 + t56 * t60;
t66 = t41 * t42;
t50 = t54 ^ 2;
t65 = t50 / t56 ^ 2;
t64 = t54 * t57;
t47 = 0.1e1 / (t55 ^ 2 * t65 + 0.1e1);
t63 = t55 * t47;
t61 = t55 * t56;
t58 = t42 ^ 2 * t41 + 0.1e1;
t51 = 0.1e1 / t56;
t40 = 0.1e1 / t43;
t39 = (0.1e1 + t65) * t63;
t36 = 0.1e1 / t38;
t35 = 0.1e1 / t58;
t34 = 0.1e1 / (t50 * t67 + 0.1e1);
t33 = t58 * t35;
t1 = [t51 * t47 * t64, t39, 0, 0, 0, 0; (-t36 * t62 - (-t45 * t50 * t51 * t63 + (t47 - 0.1e1) * t54 * t44) * t54 * t67) * t34 (t56 * t36 - (-t44 * t61 + t45 * t54 + (t44 * t56 - t45 * t62) * t39) * t54 * t37) * t57 * t34, 0, 0, 0, 0; ((-t48 * t61 - t59) * t40 - (-t49 * t61 + t60) * t66) * t35 (-t40 * t48 + t49 * t66) * t35 * t64, t33, t33, 0, 0;];
Ja_rot  = t1;
