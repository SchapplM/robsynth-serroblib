% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR7
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
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR7_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobia_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:56
% EndTime: 2019-02-26 21:31:56
% DurationCPUTime: 0.05s
% Computational Cost: add. (58->16), mult. (119->35), div. (26->9), fcn. (180->9), ass. (0->26)
t36 = cos(pkin(6));
t35 = sin(pkin(6));
t40 = cos(qJ(1));
t43 = t40 * t35;
t31 = atan2(-t43, -t36);
t29 = sin(t31);
t30 = cos(t31);
t22 = -t29 * t43 - t30 * t36;
t38 = sin(qJ(1));
t48 = 0.1e1 / t22 ^ 2 * t38 ^ 2;
t37 = sin(qJ(2));
t42 = t40 * t37;
t39 = cos(qJ(2));
t44 = t38 * t39;
t27 = t36 * t44 + t42;
t25 = 0.1e1 / t27 ^ 2;
t41 = t40 * t39;
t45 = t38 * t37;
t28 = -t36 * t45 + t41;
t47 = t28 ^ 2 * t25;
t33 = t35 ^ 2;
t32 = 0.1e1 / (0.1e1 + t40 ^ 2 * t33 / t36 ^ 2);
t46 = t32 / t36;
t24 = 0.1e1 / t27;
t23 = 0.1e1 / (0.1e1 + t47);
t1 = [-t38 * t35 * t46, 0, 0, 0, 0, 0; (-0.1e1 / t22 * t43 + (t30 * t33 * t40 * t46 + (-t32 + 0.1e1) * t35 * t29) * t35 * t48) / (t33 * t48 + 0.1e1) 0, 0, 0, 0, 0; ((-t36 * t42 - t44) * t24 - (t36 * t41 - t45) * t28 * t25) * t23 (-t24 * t27 - t47) * t23, 0, 0, 0, 0;];
Ja_rot  = t1;
