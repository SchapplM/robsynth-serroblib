% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_jacobia_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:59:51
% EndTime: 2019-02-26 20:59:51
% DurationCPUTime: 0.07s
% Computational Cost: add. (87->20), mult. (227->55), div. (53->9), fcn. (337->9), ass. (0->36)
t46 = sin(qJ(1));
t62 = t46 ^ 2;
t45 = sin(qJ(3));
t48 = cos(qJ(3));
t49 = cos(qJ(1));
t50 = t49 * t48;
t39 = atan2(t50, -t45);
t37 = sin(t39);
t38 = cos(t39);
t29 = t37 * t50 - t38 * t45;
t28 = 0.1e1 / t29 ^ 2;
t61 = t28 * t48;
t44 = sin(qJ(4));
t52 = t49 * t44;
t47 = cos(qJ(4));
t55 = t46 * t47;
t36 = t45 * t55 + t52;
t34 = 0.1e1 / t36 ^ 2;
t51 = t49 * t47;
t56 = t46 * t44;
t35 = -t45 * t56 + t51;
t60 = t35 ^ 2 * t34;
t59 = t34 * t35;
t58 = t37 * t45;
t43 = t48 ^ 2;
t57 = 0.1e1 / t45 ^ 2 * t43;
t54 = t46 * t48;
t40 = 0.1e1 / (t49 ^ 2 * t57 + 0.1e1);
t53 = t49 * t40;
t41 = 0.1e1 / t45;
t33 = 0.1e1 / t36;
t31 = (0.1e1 + t57) * t53;
t30 = 0.1e1 / (0.1e1 + t60);
t27 = 0.1e1 / t29;
t26 = 0.1e1 / (t62 * t43 * t28 + 0.1e1);
t1 = [t41 * t40 * t54, 0, t31, 0, 0, 0; (t27 * t50 - (t38 * t41 * t43 * t53 + (t40 - 0.1e1) * t48 * t37) * t62 * t61) * t26, 0 (-t45 * t27 - (-t49 * t58 - t38 * t48 + (t38 * t50 + t58) * t31) * t61) * t46 * t26, 0, 0, 0; ((-t45 * t52 - t55) * t33 - (t45 * t51 - t56) * t59) * t30, 0 (-t33 * t44 - t47 * t59) * t30 * t54 (-t33 * t36 - t60) * t30, 0, 0;];
Ja_rot  = t1;
