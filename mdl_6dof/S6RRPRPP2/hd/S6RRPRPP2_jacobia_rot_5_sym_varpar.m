% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_rot = S6RRPRPP2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:35:21
% EndTime: 2019-02-26 21:35:21
% DurationCPUTime: 0.11s
% Computational Cost: add. (353->27), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->41)
t56 = qJ(2) + pkin(9);
t54 = sin(t56);
t76 = t54 ^ 2;
t55 = cos(t56);
t62 = cos(qJ(4));
t63 = cos(qJ(1));
t65 = t63 * t62;
t60 = sin(qJ(4));
t61 = sin(qJ(1));
t68 = t61 * t60;
t44 = t55 * t68 + t65;
t70 = t54 * t60;
t41 = atan2(-t44, t70);
t37 = sin(t41);
t38 = cos(t41);
t36 = -t37 * t44 + t38 * t70;
t35 = 0.1e1 / t36 ^ 2;
t66 = t63 * t60;
t67 = t61 * t62;
t47 = t55 * t66 - t67;
t75 = t35 * t47;
t73 = t38 * t44;
t72 = t47 ^ 2 * t35;
t52 = 0.1e1 / t54;
t57 = 0.1e1 / t60;
t71 = t52 * t57;
t69 = t54 * t63;
t48 = t55 * t65 + t68;
t43 = 0.1e1 / t48 ^ 2;
t64 = t63 ^ 2 * t76 * t43;
t58 = 0.1e1 / t60 ^ 2;
t53 = 0.1e1 / t76;
t46 = t55 * t67 - t66;
t42 = 0.1e1 / t48;
t40 = 0.1e1 / (t44 ^ 2 * t53 * t58 + 0.1e1);
t39 = 0.1e1 / (0.1e1 + t64);
t34 = 0.1e1 / t36;
t33 = (t44 * t53 * t55 * t57 + t61) * t40;
t32 = 0.1e1 / (0.1e1 + t72);
t31 = (t44 * t58 * t62 - t46 * t57) * t52 * t40;
t1 = [-t47 * t40 * t71, t33, 0, t31, 0, 0; (-t44 * t34 - (-t37 + (t71 * t73 + t37) * t40) * t72) * t32 (t33 * t73 * t75 + (-t34 * t69 - (t38 * t55 + (-t33 + t61) * t54 * t37) * t75) * t60) * t32, 0 (t48 * t34 - (t38 * t54 * t62 - t37 * t46 + (-t37 * t70 - t73) * t31) * t75) * t32, 0, 0; (-t43 * t46 * t63 + t42 * t61) * t54 * t39 (-t42 * t55 * t63 - t62 * t64) * t39, 0, -t47 * t43 * t39 * t69, 0, 0;];
Ja_rot  = t1;
