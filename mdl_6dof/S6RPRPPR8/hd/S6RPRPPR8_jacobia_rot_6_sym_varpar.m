% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_jacobia_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:43:07
% EndTime: 2019-02-26 20:43:07
% DurationCPUTime: 0.07s
% Computational Cost: add. (64->19), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->34)
t43 = cos(qJ(3));
t40 = sin(qJ(3));
t44 = cos(qJ(1));
t47 = t44 * t40;
t34 = atan2(t47, t43);
t31 = sin(t34);
t32 = cos(t34);
t23 = t31 * t47 + t32 * t43;
t22 = 0.1e1 / t23 ^ 2;
t41 = sin(qJ(1));
t55 = t22 * t41 ^ 2;
t42 = cos(qJ(6));
t39 = sin(qJ(6));
t48 = t44 * t39;
t50 = t41 * t43;
t30 = -t42 * t50 - t48;
t27 = 0.1e1 / t30 ^ 2;
t46 = t44 * t42;
t28 = t39 * t50 - t46;
t54 = t27 * t28;
t53 = t31 * t43;
t36 = t40 ^ 2;
t52 = t36 / t43 ^ 2;
t51 = t40 * t41;
t33 = 0.1e1 / (t44 ^ 2 * t52 + 0.1e1);
t49 = t44 * t33;
t45 = t28 ^ 2 * t27 + 0.1e1;
t37 = 0.1e1 / t43;
t26 = 0.1e1 / t30;
t25 = (0.1e1 + t52) * t49;
t24 = 0.1e1 / t45;
t21 = 0.1e1 / t23;
t20 = 0.1e1 / (t36 * t55 + 0.1e1);
t1 = [-t37 * t33 * t51, 0, t25, 0, 0, 0; (t21 * t47 - (-t32 * t36 * t37 * t49 + (t33 - 0.1e1) * t40 * t31) * t40 * t55) * t20, 0 (t43 * t21 - (t44 * t53 - t32 * t40 + (t32 * t47 - t53) * t25) * t40 * t22) * t41 * t20, 0, 0, 0; ((-t41 * t42 - t43 * t48) * t26 + (t41 * t39 - t43 * t46) * t54) * t24, 0 (t26 * t39 + t42 * t54) * t24 * t51, 0, 0, t45 * t24;];
Ja_rot  = t1;
