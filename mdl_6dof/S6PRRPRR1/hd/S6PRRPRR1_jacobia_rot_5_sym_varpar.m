% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:04:06
% EndTime: 2019-02-26 20:04:06
% DurationCPUTime: 0.09s
% Computational Cost: add. (206->19), mult. (316->43), div. (52->11), fcn. (458->11), ass. (0->31)
t62 = sin(pkin(11));
t63 = sin(pkin(6));
t72 = t62 * t63;
t67 = cos(qJ(2));
t71 = t63 * t67;
t65 = cos(pkin(6));
t66 = sin(qJ(2));
t70 = t65 * t66;
t69 = t65 * t67;
t64 = cos(pkin(11));
t55 = -t62 * t70 + t64 * t67;
t59 = qJ(3) + pkin(12) + qJ(5);
t57 = sin(t59);
t58 = cos(t59);
t46 = t55 * t58 + t57 * t72;
t44 = 0.1e1 / t46 ^ 2;
t45 = t55 * t57 - t58 * t72;
t68 = t45 ^ 2 * t44 + 0.1e1;
t61 = 0.1e1 / t67 ^ 2;
t54 = t62 * t69 + t64 * t66;
t53 = t62 * t67 + t64 * t70;
t51 = t62 * t66 - t64 * t69;
t49 = atan2(-t51, -t71);
t48 = cos(t49);
t47 = sin(t49);
t43 = 0.1e1 / t68;
t42 = -t47 * t51 - t48 * t71;
t41 = 0.1e1 / t42 ^ 2;
t39 = (t53 / t67 + t66 * t51 * t61) / t63 / (0.1e1 + t51 ^ 2 / t63 ^ 2 * t61);
t38 = t68 * t43;
t1 = [0, t39, 0, 0, 0, 0; 0 (t55 / t42 - (t48 * t63 * t66 - t47 * t53 + (t47 * t71 - t48 * t51) * t39) * t54 * t41) / (t54 ^ 2 * t41 + 0.1e1) 0, 0, 0, 0; 0 (-t57 / t46 + t58 * t45 * t44) * t54 * t43, t38, 0, t38, 0;];
Ja_rot  = t1;
