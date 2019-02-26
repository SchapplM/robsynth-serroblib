% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:31:47
% EndTime: 2019-02-26 22:31:47
% DurationCPUTime: 0.06s
% Computational Cost: add. (632->19), mult. (313->39), div. (91->9), fcn. (471->7), ass. (0->30)
t62 = cos(qJ(1));
t60 = t62 ^ 2;
t57 = qJ(2) + qJ(3) + qJ(4);
t56 = cos(t57);
t55 = sin(t57);
t61 = sin(qJ(1));
t66 = t61 * t55;
t49 = atan2(-t66, -t56);
t47 = sin(t49);
t48 = cos(t49);
t45 = -t47 * t66 - t48 * t56;
t44 = 0.1e1 / t45 ^ 2;
t72 = t44 * t55;
t71 = t47 * t56;
t52 = t55 ^ 2;
t63 = t56 ^ 2;
t70 = t52 / t63;
t69 = t55 * t62;
t64 = t61 ^ 2;
t68 = 0.1e1 / t64 * t60;
t50 = 0.1e1 / (t64 * t70 + 0.1e1);
t67 = t61 * t50;
t51 = 0.1e1 / (t63 * t68 + 0.1e1);
t65 = 0.1e1 / t61 * t51 * t69;
t53 = 0.1e1 / t56;
t46 = (0.1e1 + t70) * t67;
t43 = 0.1e1 / t45;
t42 = 0.1e1 / (t60 * t52 * t44 + 0.1e1);
t41 = (t56 * t43 - (-t61 * t71 + t48 * t55 + (-t48 * t66 + t71) * t46) * t72) * t62 * t42;
t1 = [t53 * t50 * t69, t46, t46, t46, 0, 0; (-t43 * t66 - (-t48 * t52 * t53 * t67 + (t50 - 0.1e1) * t55 * t47) * t60 * t72) * t42, t41, t41, t41, 0, 0; (-0.1e1 - t68) * t56 * t51, -t65, -t65, -t65, 0, 0;];
Ja_rot  = t1;
