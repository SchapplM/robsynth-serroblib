% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR2_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobia_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:06
% EndTime: 2019-02-26 20:11:06
% DurationCPUTime: 0.09s
% Computational Cost: add. (162->19), mult. (316->43), div. (52->11), fcn. (458->11), ass. (0->31)
t59 = sin(pkin(11));
t60 = sin(pkin(6));
t69 = t59 * t60;
t64 = cos(qJ(2));
t68 = t60 * t64;
t62 = cos(pkin(6));
t63 = sin(qJ(2));
t67 = t62 * t63;
t66 = t62 * t64;
t61 = cos(pkin(11));
t52 = -t59 * t67 + t61 * t64;
t58 = qJ(3) + qJ(4);
t54 = sin(t58);
t55 = cos(t58);
t43 = t52 * t55 + t54 * t69;
t41 = 0.1e1 / t43 ^ 2;
t42 = t52 * t54 - t55 * t69;
t65 = t42 ^ 2 * t41 + 0.1e1;
t57 = 0.1e1 / t64 ^ 2;
t51 = t59 * t66 + t61 * t63;
t50 = t59 * t64 + t61 * t67;
t48 = t59 * t63 - t61 * t66;
t46 = atan2(-t48, -t68);
t45 = cos(t46);
t44 = sin(t46);
t40 = 0.1e1 / t65;
t39 = -t44 * t48 - t45 * t68;
t38 = 0.1e1 / t39 ^ 2;
t36 = (t50 / t64 + t63 * t48 * t57) / t60 / (0.1e1 + t48 ^ 2 / t60 ^ 2 * t57);
t35 = t65 * t40;
t1 = [0, t36, 0, 0, 0, 0; 0 (t52 / t39 - (t45 * t60 * t63 - t44 * t50 + (t44 * t68 - t45 * t48) * t36) * t51 * t38) / (t51 ^ 2 * t38 + 0.1e1) 0, 0, 0, 0; 0 (-t54 / t43 + t55 * t42 * t41) * t51 * t40, t35, t35, 0, 0;];
Ja_rot  = t1;
