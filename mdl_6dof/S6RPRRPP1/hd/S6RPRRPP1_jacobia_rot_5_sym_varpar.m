% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:56:16
% EndTime: 2019-02-26 20:56:17
% DurationCPUTime: 0.11s
% Computational Cost: add. (266->22), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
t47 = qJ(1) + pkin(9);
t45 = cos(t47);
t62 = t45 ^ 2;
t52 = cos(qJ(3));
t43 = sin(t47);
t51 = sin(qJ(3));
t58 = t43 * t51;
t40 = atan2(-t58, -t52);
t38 = sin(t40);
t39 = cos(t40);
t32 = -t38 * t58 - t39 * t52;
t31 = 0.1e1 / t32 ^ 2;
t61 = t31 * t51;
t46 = qJ(4) + pkin(10);
t42 = sin(t46);
t44 = cos(t46);
t55 = t45 * t52;
t37 = t43 * t42 + t44 * t55;
t35 = 0.1e1 / t37 ^ 2;
t36 = t42 * t55 - t43 * t44;
t60 = t35 * t36;
t48 = t51 ^ 2;
t54 = t48 / t52 ^ 2;
t41 = 0.1e1 / (t43 ^ 2 * t54 + 0.1e1);
t59 = t43 * t41;
t57 = t43 * t52;
t56 = t45 * t51;
t53 = t36 ^ 2 * t35 + 0.1e1;
t49 = 0.1e1 / t52;
t34 = 0.1e1 / t37;
t33 = (0.1e1 + t54) * t59;
t30 = 0.1e1 / t32;
t29 = 0.1e1 / t53;
t28 = 0.1e1 / (t62 * t48 * t31 + 0.1e1);
t1 = [t49 * t41 * t56, 0, t33, 0, 0, 0; (-t30 * t58 - (-t39 * t48 * t49 * t59 + (t41 - 0.1e1) * t51 * t38) * t62 * t61) * t28, 0 (t52 * t30 - (-t38 * t57 + t39 * t51 + (t38 * t52 - t39 * t58) * t33) * t61) * t45 * t28, 0, 0, 0; ((-t42 * t57 - t45 * t44) * t34 - (t45 * t42 - t44 * t57) * t60) * t29, 0 (-t34 * t42 + t44 * t60) * t29 * t56, t53 * t29, 0, 0;];
Ja_rot  = t1;
