% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR5_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobia_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:47
% EndTime: 2019-02-26 21:30:48
% DurationCPUTime: 0.05s
% Computational Cost: add. (58->16), mult. (119->35), div. (26->9), fcn. (180->9), ass. (0->26)
t35 = cos(pkin(6));
t34 = sin(pkin(6));
t39 = cos(qJ(1));
t42 = t39 * t34;
t30 = atan2(-t42, -t35);
t28 = sin(t30);
t29 = cos(t30);
t21 = -t28 * t42 - t29 * t35;
t37 = sin(qJ(1));
t47 = 0.1e1 / t21 ^ 2 * t37 ^ 2;
t38 = cos(qJ(2));
t40 = t39 * t38;
t36 = sin(qJ(2));
t44 = t37 * t36;
t27 = -t35 * t44 + t40;
t25 = 0.1e1 / t27 ^ 2;
t41 = t39 * t36;
t43 = t37 * t38;
t26 = -t35 * t43 - t41;
t46 = t26 ^ 2 * t25;
t32 = t34 ^ 2;
t31 = 0.1e1 / (0.1e1 + t39 ^ 2 * t32 / t35 ^ 2);
t45 = t31 / t35;
t24 = 0.1e1 / t27;
t22 = 0.1e1 / (0.1e1 + t46);
t1 = [-t37 * t34 * t45, 0, 0, 0, 0, 0; (-0.1e1 / t21 * t42 + (t29 * t32 * t39 * t45 + (-t31 + 0.1e1) * t34 * t28) * t34 * t47) / (t32 * t47 + 0.1e1) 0, 0, 0, 0, 0; ((-t35 * t40 + t44) * t24 - (-t35 * t41 - t43) * t26 * t25) * t22 (-t24 * t27 - t46) * t22, 0, 0, 0, 0;];
Ja_rot  = t1;
