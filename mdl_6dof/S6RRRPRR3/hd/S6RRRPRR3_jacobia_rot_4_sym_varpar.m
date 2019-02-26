% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR3_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_jacobia_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:17:07
% EndTime: 2019-02-26 22:17:07
% DurationCPUTime: 0.06s
% Computational Cost: add. (299->18), mult. (227->41), div. (75->9), fcn. (351->7), ass. (0->29)
t54 = cos(qJ(1));
t55 = t54 ^ 2;
t52 = qJ(2) + qJ(3);
t48 = cos(t52);
t47 = sin(t52);
t53 = sin(qJ(1));
t57 = t53 * t47;
t41 = atan2(-t57, -t48);
t39 = sin(t41);
t40 = cos(t41);
t37 = -t39 * t57 - t40 * t48;
t36 = 0.1e1 / t37 ^ 2;
t62 = t36 * t47;
t61 = t39 * t48;
t44 = t47 ^ 2;
t46 = 0.1e1 / t48 ^ 2;
t60 = t44 * t46;
t49 = t53 ^ 2;
t59 = t49 / t55;
t42 = 0.1e1 / (t49 * t60 + 0.1e1);
t58 = t53 * t42;
t43 = 0.1e1 / (t46 * t59 + 0.1e1);
t56 = 0.1e1 / t54 * t46 * t43 * t57;
t45 = 0.1e1 / t48;
t38 = (0.1e1 + t60) * t58;
t35 = 0.1e1 / t37;
t34 = 0.1e1 / (t55 * t44 * t36 + 0.1e1);
t33 = (t48 * t35 - (-t53 * t61 + t40 * t47 + (-t40 * t57 + t61) * t38) * t62) * t54 * t34;
t1 = [t54 * t47 * t45 * t42, t38, t38, 0, 0, 0; (-t35 * t57 - (-t40 * t44 * t45 * t58 + (t42 - 0.1e1) * t47 * t39) * t55 * t62) * t34, t33, t33, 0, 0, 0; (-0.1e1 - t59) * t45 * t43, -t56, -t56, 0, 0, 0;];
Ja_rot  = t1;
