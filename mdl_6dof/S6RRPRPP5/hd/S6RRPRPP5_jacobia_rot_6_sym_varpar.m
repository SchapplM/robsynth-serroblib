% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPP5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_jacobia_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:06
% EndTime: 2019-02-26 21:37:06
% DurationCPUTime: 0.07s
% Computational Cost: add. (88->21), mult. (227->54), div. (53->9), fcn. (337->9), ass. (0->35)
t44 = sin(qJ(2));
t45 = sin(qJ(1));
t47 = cos(qJ(2));
t52 = t45 * t47;
t38 = atan2(t52, -t44);
t36 = sin(t38);
t37 = cos(t38);
t28 = t36 * t52 - t37 * t44;
t27 = 0.1e1 / t28 ^ 2;
t48 = cos(qJ(1));
t60 = t27 * t48 ^ 2;
t43 = sin(qJ(4));
t50 = t48 * t43;
t46 = cos(qJ(4));
t53 = t45 * t46;
t35 = t44 * t50 + t53;
t33 = 0.1e1 / t35 ^ 2;
t49 = t48 * t46;
t54 = t45 * t43;
t34 = t44 * t49 - t54;
t59 = t34 ^ 2 * t33;
t58 = t33 * t34;
t57 = t36 * t44;
t42 = t47 ^ 2;
t56 = 0.1e1 / t44 ^ 2 * t42;
t39 = 0.1e1 / (t45 ^ 2 * t56 + 0.1e1);
t55 = t45 * t39;
t51 = t47 * t48;
t40 = 0.1e1 / t44;
t32 = 0.1e1 / t35;
t30 = (0.1e1 + t56) * t55;
t29 = 0.1e1 / (0.1e1 + t59);
t26 = 0.1e1 / t28;
t25 = 0.1e1 / (t42 * t60 + 0.1e1);
t1 = [-t40 * t39 * t51, t30, 0, 0, 0, 0; (t26 * t52 + (-t37 * t40 * t42 * t55 + (-t39 + 0.1e1) * t47 * t36) * t47 * t60) * t25 (t44 * t26 + (-t45 * t57 - t37 * t47 + (t37 * t52 + t57) * t30) * t47 * t27) * t48 * t25, 0, 0, 0, 0; ((-t44 * t53 - t50) * t32 - (-t44 * t54 + t49) * t58) * t29 (t32 * t46 - t43 * t58) * t29 * t51, 0 (-t32 * t35 - t59) * t29, 0, 0;];
Ja_rot  = t1;
