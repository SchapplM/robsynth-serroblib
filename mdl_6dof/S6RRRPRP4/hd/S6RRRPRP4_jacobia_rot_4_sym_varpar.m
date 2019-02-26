% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP4_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_jacobia_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:11:05
% EndTime: 2019-02-26 22:11:05
% DurationCPUTime: 0.06s
% Computational Cost: add. (297->18), mult. (231->39), div. (67->9), fcn. (349->7), ass. (0->30)
t55 = cos(qJ(1));
t52 = t55 ^ 2;
t53 = qJ(2) + qJ(3);
t49 = cos(t53);
t48 = sin(t53);
t54 = sin(qJ(1));
t59 = t54 * t48;
t42 = atan2(-t59, -t49);
t40 = sin(t42);
t41 = cos(t42);
t38 = -t40 * t59 - t41 * t49;
t37 = 0.1e1 / t38 ^ 2;
t65 = t37 * t48;
t64 = t40 * t49;
t45 = t48 ^ 2;
t56 = t49 ^ 2;
t63 = t45 / t56;
t62 = t48 * t55;
t57 = t54 ^ 2;
t61 = 0.1e1 / t57 * t52;
t43 = 0.1e1 / (t57 * t63 + 0.1e1);
t60 = t54 * t43;
t44 = 0.1e1 / (t56 * t61 + 0.1e1);
t58 = 0.1e1 / t54 * t44 * t62;
t46 = 0.1e1 / t49;
t39 = (0.1e1 + t63) * t60;
t36 = 0.1e1 / t38;
t35 = 0.1e1 / (t52 * t45 * t37 + 0.1e1);
t34 = (t49 * t36 - (-t54 * t64 + t41 * t48 + (-t41 * t59 + t64) * t39) * t65) * t55 * t35;
t1 = [t46 * t43 * t62, t39, t39, 0, 0, 0; (-t36 * t59 - (-t41 * t45 * t46 * t60 + (t43 - 0.1e1) * t48 * t40) * t52 * t65) * t35, t34, t34, 0, 0, 0; (-0.1e1 - t61) * t49 * t44, -t58, -t58, 0, 0, 0;];
Ja_rot  = t1;
