% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_jacobia_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:28:35
% EndTime: 2019-02-26 20:28:35
% DurationCPUTime: 0.07s
% Computational Cost: add. (87->20), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->33)
t43 = cos(qJ(4));
t40 = sin(qJ(4));
t41 = sin(qJ(1));
t49 = t41 * t40;
t33 = atan2(-t49, t43);
t31 = sin(t33);
t32 = cos(t33);
t23 = -t31 * t49 + t32 * t43;
t22 = 0.1e1 / t23 ^ 2;
t44 = cos(qJ(1));
t54 = t22 * t44 ^ 2;
t42 = cos(qJ(6));
t39 = sin(qJ(6));
t47 = t44 * t39;
t30 = -t41 * t42 - t43 * t47;
t27 = 0.1e1 / t30 ^ 2;
t46 = t44 * t42;
t28 = t41 * t39 - t43 * t46;
t53 = t27 * t28;
t36 = t40 ^ 2;
t52 = t36 / t43 ^ 2;
t51 = t40 * t44;
t34 = 0.1e1 / (t41 ^ 2 * t52 + 0.1e1);
t50 = t41 * t34;
t48 = t41 * t43;
t45 = t28 ^ 2 * t27 + 0.1e1;
t37 = 0.1e1 / t43;
t26 = 0.1e1 / t30;
t25 = (-0.1e1 - t52) * t50;
t24 = 0.1e1 / t45;
t21 = 0.1e1 / t23;
t20 = 0.1e1 / (t36 * t54 + 0.1e1);
t1 = [-t37 * t34 * t51, 0, 0, t25, 0, 0; (-t21 * t49 - (t32 * t36 * t37 * t50 + (t34 - 0.1e1) * t40 * t31) * t40 * t54) * t20, 0, 0 (t43 * t21 - (-t31 * t48 - t32 * t40 + (-t31 * t43 - t32 * t49) * t25) * t40 * t22) * t44 * t20, 0, 0; ((-t42 * t48 - t47) * t26 + (t39 * t48 - t46) * t53) * t24, 0, 0 (-t26 * t42 + t39 * t53) * t24 * t51, 0, t45 * t24;];
Ja_rot  = t1;
