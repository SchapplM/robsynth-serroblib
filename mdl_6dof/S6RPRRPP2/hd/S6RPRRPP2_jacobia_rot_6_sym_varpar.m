% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:56:53
% EndTime: 2019-02-26 20:56:53
% DurationCPUTime: 0.07s
% Computational Cost: add. (175->20), mult. (227->57), div. (53->9), fcn. (337->9), ass. (0->35)
t41 = qJ(1) + pkin(9);
t40 = cos(t41);
t59 = t40 ^ 2;
t48 = cos(qJ(3));
t39 = sin(t41);
t46 = sin(qJ(3));
t53 = t39 * t46;
t38 = atan2(t53, t48);
t35 = sin(t38);
t36 = cos(t38);
t29 = t35 * t53 + t36 * t48;
t28 = 0.1e1 / t29 ^ 2;
t58 = t28 * t46;
t45 = sin(qJ(4));
t47 = cos(qJ(4));
t49 = t47 * t48;
t34 = t39 * t45 + t40 * t49;
t32 = 0.1e1 / t34 ^ 2;
t50 = t45 * t48;
t33 = t39 * t47 - t40 * t50;
t57 = t33 ^ 2 * t32;
t56 = t32 * t33;
t55 = t35 * t48;
t42 = t46 ^ 2;
t51 = t42 / t48 ^ 2;
t37 = 0.1e1 / (t39 ^ 2 * t51 + 0.1e1);
t54 = t39 * t37;
t52 = t40 * t46;
t43 = 0.1e1 / t48;
t31 = 0.1e1 / t34;
t27 = 0.1e1 / t29;
t26 = (0.1e1 + t51) * t54;
t25 = 0.1e1 / (0.1e1 + t57);
t24 = 0.1e1 / (t59 * t42 * t28 + 0.1e1);
t1 = [t43 * t37 * t52, 0, t26, 0, 0, 0; (t27 * t53 + (t36 * t42 * t43 * t54 + (-t37 + 0.1e1) * t46 * t35) * t59 * t58) * t24, 0 (-t48 * t27 + (t39 * t55 - t36 * t46 + (t36 * t53 - t55) * t26) * t58) * t40 * t24, 0, 0, 0; ((t39 * t50 + t40 * t47) * t31 - (-t39 * t49 + t40 * t45) * t56) * t25, 0 (t31 * t45 + t47 * t56) * t25 * t52 (-t31 * t34 - t57) * t25, 0, 0;];
Ja_rot  = t1;
