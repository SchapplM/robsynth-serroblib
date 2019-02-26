% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S6RRRRRP11
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRP11_jacobia_rot_2_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_jacobia_rot_2_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_jacobia_rot_2_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:38
% EndTime: 2019-02-26 22:45:38
% DurationCPUTime: 0.05s
% Computational Cost: add. (35->13), mult. (116->34), div. (25->9), fcn. (175->9), ass. (0->25)
t30 = cos(pkin(6));
t29 = sin(pkin(6));
t34 = cos(qJ(1));
t38 = t34 * t29;
t26 = atan2(t38, t30);
t23 = sin(t26);
t24 = cos(t26);
t18 = t23 * t38 + t24 * t30;
t32 = sin(qJ(1));
t42 = 0.1e1 / t18 ^ 2 * t32 ^ 2;
t27 = t29 ^ 2;
t25 = 0.1e1 / (0.1e1 + t34 ^ 2 * t27 / t30 ^ 2);
t41 = t25 / t30;
t31 = sin(qJ(2));
t40 = t32 * t31;
t33 = cos(qJ(2));
t39 = t32 * t33;
t37 = t34 * t31;
t36 = t34 * t33;
t22 = -t30 * t40 + t36;
t20 = 0.1e1 / t22 ^ 2;
t21 = t30 * t39 + t37;
t35 = t21 ^ 2 * t20 + 0.1e1;
t19 = 0.1e1 / t35;
t1 = [-t32 * t29 * t41, 0, 0, 0, 0, 0; (0.1e1 / t18 * t38 - (-t24 * t27 * t34 * t41 + (t25 - 0.1e1) * t29 * t23) * t29 * t42) / (t27 * t42 + 0.1e1) 0, 0, 0, 0, 0; ((t30 * t36 - t40) / t22 - (-t30 * t37 - t39) * t21 * t20) * t19, t35 * t19, 0, 0, 0, 0;];
Ja_rot  = t1;
