% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:38:16
% EndTime: 2019-02-26 21:38:16
% DurationCPUTime: 0.06s
% Computational Cost: add. (467->18), mult. (231->39), div. (67->9), fcn. (349->7), ass. (0->30)
t58 = cos(qJ(1));
t56 = t58 ^ 2;
t53 = qJ(2) + pkin(10) + qJ(4);
t52 = cos(t53);
t51 = sin(t53);
t57 = sin(qJ(1));
t62 = t57 * t51;
t45 = atan2(-t62, -t52);
t43 = sin(t45);
t44 = cos(t45);
t41 = -t43 * t62 - t44 * t52;
t40 = 0.1e1 / t41 ^ 2;
t68 = t40 * t51;
t67 = t43 * t52;
t48 = t51 ^ 2;
t59 = t52 ^ 2;
t66 = t48 / t59;
t65 = t51 * t58;
t60 = t57 ^ 2;
t64 = 0.1e1 / t60 * t56;
t46 = 0.1e1 / (t60 * t66 + 0.1e1);
t63 = t57 * t46;
t47 = 0.1e1 / (t59 * t64 + 0.1e1);
t61 = 0.1e1 / t57 * t47 * t65;
t49 = 0.1e1 / t52;
t42 = (0.1e1 + t66) * t63;
t39 = 0.1e1 / t41;
t38 = 0.1e1 / (t56 * t48 * t40 + 0.1e1);
t37 = (t52 * t39 - (-t57 * t67 + t44 * t51 + (-t44 * t62 + t67) * t42) * t68) * t58 * t38;
t1 = [t49 * t46 * t65, t42, 0, t42, 0, 0; (-t39 * t62 - (-t44 * t48 * t49 * t63 + (t46 - 0.1e1) * t51 * t43) * t56 * t68) * t38, t37, 0, t37, 0, 0; (-0.1e1 - t64) * t52 * t47, -t61, 0, -t61, 0, 0;];
Ja_rot  = t1;
