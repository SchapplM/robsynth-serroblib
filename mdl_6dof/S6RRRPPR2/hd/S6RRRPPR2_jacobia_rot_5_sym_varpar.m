% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:03:56
% EndTime: 2019-02-26 22:03:56
% DurationCPUTime: 0.06s
% Computational Cost: add. (467->18), mult. (231->39), div. (67->9), fcn. (349->7), ass. (0->30)
t61 = cos(qJ(1));
t59 = t61 ^ 2;
t56 = qJ(2) + qJ(3) + pkin(10);
t55 = cos(t56);
t54 = sin(t56);
t60 = sin(qJ(1));
t65 = t60 * t54;
t48 = atan2(-t65, -t55);
t46 = sin(t48);
t47 = cos(t48);
t44 = -t46 * t65 - t47 * t55;
t43 = 0.1e1 / t44 ^ 2;
t71 = t43 * t54;
t70 = t46 * t55;
t51 = t54 ^ 2;
t62 = t55 ^ 2;
t69 = t51 / t62;
t68 = t54 * t61;
t63 = t60 ^ 2;
t67 = 0.1e1 / t63 * t59;
t49 = 0.1e1 / (t63 * t69 + 0.1e1);
t66 = t60 * t49;
t50 = 0.1e1 / (t62 * t67 + 0.1e1);
t64 = 0.1e1 / t60 * t50 * t68;
t52 = 0.1e1 / t55;
t45 = (0.1e1 + t69) * t66;
t42 = 0.1e1 / t44;
t41 = 0.1e1 / (t59 * t51 * t43 + 0.1e1);
t40 = (t55 * t42 - (-t60 * t70 + t47 * t54 + (-t47 * t65 + t70) * t45) * t71) * t61 * t41;
t1 = [t52 * t49 * t68, t45, t45, 0, 0, 0; (-t42 * t65 - (-t47 * t51 * t52 * t66 + (t49 - 0.1e1) * t54 * t46) * t59 * t71) * t41, t40, t40, 0, 0, 0; (-0.1e1 - t67) * t55 * t50, -t64, -t64, 0, 0, 0;];
Ja_rot  = t1;
