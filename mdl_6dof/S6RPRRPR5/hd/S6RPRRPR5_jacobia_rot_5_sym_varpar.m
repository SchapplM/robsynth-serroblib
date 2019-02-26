% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:03:12
% EndTime: 2019-02-26 21:03:12
% DurationCPUTime: 0.06s
% Computational Cost: add. (467->18), mult. (231->39), div. (67->9), fcn. (349->7), ass. (0->30)
t56 = cos(qJ(1));
t54 = t56 ^ 2;
t51 = pkin(10) + qJ(3) + qJ(4);
t50 = cos(t51);
t49 = sin(t51);
t55 = sin(qJ(1));
t60 = t55 * t49;
t43 = atan2(-t60, -t50);
t41 = sin(t43);
t42 = cos(t43);
t39 = -t41 * t60 - t42 * t50;
t38 = 0.1e1 / t39 ^ 2;
t66 = t38 * t49;
t65 = t41 * t50;
t46 = t49 ^ 2;
t57 = t50 ^ 2;
t64 = t46 / t57;
t63 = t49 * t56;
t58 = t55 ^ 2;
t62 = 0.1e1 / t58 * t54;
t44 = 0.1e1 / (t58 * t64 + 0.1e1);
t61 = t55 * t44;
t45 = 0.1e1 / (t57 * t62 + 0.1e1);
t59 = 0.1e1 / t55 * t45 * t63;
t47 = 0.1e1 / t50;
t40 = (0.1e1 + t64) * t61;
t37 = 0.1e1 / t39;
t36 = 0.1e1 / (t54 * t46 * t38 + 0.1e1);
t35 = (t50 * t37 - (-t55 * t65 + t42 * t49 + (-t42 * t60 + t65) * t40) * t66) * t56 * t36;
t1 = [t47 * t44 * t63, 0, t40, t40, 0, 0; (-t37 * t60 - (-t42 * t46 * t47 * t61 + (t44 - 0.1e1) * t49 * t41) * t54 * t66) * t36, 0, t35, t35, 0, 0; (-0.1e1 - t62) * t50 * t45, 0, -t59, -t59, 0, 0;];
Ja_rot  = t1;
