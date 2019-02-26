% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR6_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:41:58
% EndTime: 2019-02-26 20:41:58
% DurationCPUTime: 0.07s
% Computational Cost: add. (198->20), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->35)
t38 = qJ(3) + pkin(9);
t36 = sin(t38);
t37 = cos(t38);
t42 = cos(qJ(1));
t45 = t42 * t37;
t31 = atan2(-t45, t36);
t29 = sin(t31);
t30 = cos(t31);
t22 = -t29 * t45 + t30 * t36;
t21 = 0.1e1 / t22 ^ 2;
t41 = sin(qJ(1));
t53 = t21 * t41 ^ 2;
t39 = sin(pkin(10));
t44 = t42 * t39;
t40 = cos(pkin(10));
t47 = t41 * t40;
t28 = t36 * t47 + t44;
t26 = 0.1e1 / t28 ^ 2;
t43 = t42 * t40;
t48 = t41 * t39;
t27 = t36 * t48 - t43;
t52 = t26 * t27;
t51 = t29 * t36;
t35 = t37 ^ 2;
t50 = 0.1e1 / t36 ^ 2 * t35;
t49 = t37 * t41;
t32 = 0.1e1 / (t42 ^ 2 * t50 + 0.1e1);
t46 = t42 * t32;
t33 = 0.1e1 / t36;
t25 = 0.1e1 / t28;
t24 = 0.1e1 / (t27 ^ 2 * t26 + 0.1e1);
t23 = (0.1e1 + t50) * t46;
t20 = 0.1e1 / t22;
t19 = 0.1e1 / (t35 * t53 + 0.1e1);
t1 = [t33 * t32 * t49, 0, t23, 0, 0, 0; (-t20 * t45 + (-t30 * t33 * t35 * t46 + (-t32 + 0.1e1) * t37 * t29) * t37 * t53) * t19, 0 (t36 * t20 + (t42 * t51 + t30 * t37 + (-t30 * t45 - t51) * t23) * t37 * t21) * t41 * t19, 0, 0, 0; ((t36 * t44 + t47) * t25 - (t36 * t43 - t48) * t52) * t24, 0 (t25 * t39 - t40 * t52) * t24 * t49, 0, 0, 0;];
Ja_rot  = t1;
