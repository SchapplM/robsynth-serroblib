% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP11_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:51:38
% EndTime: 2019-02-26 21:51:38
% DurationCPUTime: 0.10s
% Computational Cost: add. (159->21), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->37)
t56 = sin(qJ(2));
t57 = sin(qJ(1));
t58 = cos(qJ(2));
t64 = t57 * t58;
t48 = atan2(-t64, t56);
t46 = sin(t48);
t47 = cos(t48);
t40 = -t46 * t64 + t47 * t56;
t39 = 0.1e1 / t40 ^ 2;
t59 = cos(qJ(1));
t71 = t39 * t59 ^ 2;
t55 = qJ(4) + qJ(5);
t50 = sin(t55);
t62 = t59 * t50;
t51 = cos(t55);
t65 = t57 * t51;
t45 = t56 * t62 + t65;
t43 = 0.1e1 / t45 ^ 2;
t61 = t59 * t51;
t66 = t57 * t50;
t44 = -t56 * t61 + t66;
t70 = t43 * t44;
t69 = t46 * t56;
t54 = t58 ^ 2;
t68 = 0.1e1 / t56 ^ 2 * t54;
t49 = 0.1e1 / (t57 ^ 2 * t68 + 0.1e1);
t67 = t57 * t49;
t63 = t58 * t59;
t60 = t44 ^ 2 * t43 + 0.1e1;
t52 = 0.1e1 / t56;
t42 = 0.1e1 / t45;
t41 = (0.1e1 + t68) * t67;
t38 = 0.1e1 / t40;
t37 = 0.1e1 / t60;
t36 = 0.1e1 / (t54 * t71 + 0.1e1);
t35 = t60 * t37;
t1 = [-t52 * t49 * t63, t41, 0, 0, 0, 0; (-t38 * t64 - (t47 * t52 * t54 * t67 + (t49 - 0.1e1) * t58 * t46) * t58 * t71) * t36 (-t56 * t38 - (t57 * t69 + t47 * t58 + (-t47 * t64 - t69) * t41) * t58 * t39) * t59 * t36, 0, 0, 0, 0; ((t56 * t65 + t62) * t42 - (-t56 * t66 + t61) * t70) * t37 (-t42 * t51 - t50 * t70) * t37 * t63, 0, t35, t35, 0;];
Ja_rot  = t1;
