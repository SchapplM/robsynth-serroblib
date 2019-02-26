% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPPRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPPRR2_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:22
% EndTime: 2019-02-26 19:45:22
% DurationCPUTime: 0.06s
% Computational Cost: add. (161->15), mult. (461->36), div. (29->11), fcn. (658->11), ass. (0->26)
t49 = cos(pkin(6));
t44 = sin(pkin(11));
t47 = cos(pkin(11));
t50 = sin(qJ(2));
t51 = cos(qJ(2));
t53 = t51 * t44 + t50 * t47;
t41 = t53 * t49;
t42 = t50 * t44 - t51 * t47;
t45 = sin(pkin(10));
t48 = cos(pkin(10));
t36 = -t45 * t41 - t48 * t42;
t52 = t42 * t49;
t46 = sin(pkin(6));
t40 = t53 * t46;
t39 = t42 * t46;
t38 = 0.1e1 / t39 ^ 2;
t34 = t45 * t52 - t48 * t53;
t33 = -t45 * t53 - t48 * t52;
t32 = -t48 * t41 + t45 * t42;
t31 = atan2(t33, t39);
t29 = cos(t31);
t28 = sin(t31);
t27 = t28 * t33 + t29 * t39;
t26 = 0.1e1 / t27 ^ 2;
t24 = (t32 / t39 - t40 * t33 * t38) / (t33 ^ 2 * t38 + 0.1e1);
t1 = [0, t24, 0, 0, 0, 0; 0 (t36 / t27 + (t28 * t32 + t29 * t40 + (-t28 * t39 + t29 * t33) * t24) * t34 * t26) / (t34 ^ 2 * t26 + 0.1e1) 0, 0, 0, 0; 0, t34 / t45 / t46 / (0.1e1 + t36 ^ 2 / t45 ^ 2 / t46 ^ 2) 0, 0, 0, 0;];
Ja_rot  = t1;
