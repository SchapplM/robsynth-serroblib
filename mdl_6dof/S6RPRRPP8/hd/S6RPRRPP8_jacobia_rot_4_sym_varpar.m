% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRPP8
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
% Datum: 2019-02-26 21:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP8_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_jacobia_rot_4_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:00:25
% EndTime: 2019-02-26 21:00:25
% DurationCPUTime: 0.06s
% Computational Cost: add. (86->19), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->36)
t40 = sin(qJ(1));
t56 = t40 ^ 2;
t39 = sin(qJ(3));
t42 = cos(qJ(3));
t43 = cos(qJ(1));
t45 = t43 * t42;
t33 = atan2(-t45, t39);
t31 = sin(t33);
t32 = cos(t33);
t24 = -t31 * t45 + t32 * t39;
t23 = 0.1e1 / t24 ^ 2;
t55 = t23 * t42;
t38 = sin(qJ(4));
t47 = t43 * t38;
t41 = cos(qJ(4));
t50 = t40 * t41;
t30 = t39 * t50 + t47;
t28 = 0.1e1 / t30 ^ 2;
t46 = t43 * t41;
t51 = t40 * t38;
t29 = t39 * t51 - t46;
t54 = t28 * t29;
t53 = t31 * t39;
t37 = t42 ^ 2;
t52 = 0.1e1 / t39 ^ 2 * t37;
t49 = t40 * t42;
t34 = 0.1e1 / (t43 ^ 2 * t52 + 0.1e1);
t48 = t43 * t34;
t44 = t29 ^ 2 * t28 + 0.1e1;
t35 = 0.1e1 / t39;
t27 = 0.1e1 / t30;
t26 = (0.1e1 + t52) * t48;
t25 = 0.1e1 / t44;
t22 = 0.1e1 / t24;
t21 = 0.1e1 / (t56 * t37 * t23 + 0.1e1);
t1 = [t35 * t34 * t49, 0, t26, 0, 0, 0; (-t22 * t45 + (-t32 * t35 * t37 * t48 + (-t34 + 0.1e1) * t42 * t31) * t56 * t55) * t21, 0 (t39 * t22 + (t43 * t53 + t32 * t42 + (-t32 * t45 - t53) * t26) * t55) * t40 * t21, 0, 0, 0; ((t39 * t47 + t50) * t27 - (t39 * t46 - t51) * t54) * t25, 0 (t27 * t38 - t41 * t54) * t25 * t49, t44 * t25, 0, 0;];
Ja_rot  = t1;
