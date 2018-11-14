% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S4RPPP1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
%
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S4RPPP1_jacobia_rot_3_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobia_rot_3_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobia_rot_3_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:47
% EndTime: 2018-11-14 13:45:47
% DurationCPUTime: 0.06s
% Computational Cost: add. (199->19), mult. (196->29), div. (24->11), fcn. (226->12), ass. (0->24)
t39 = sin(pkin(6));
t42 = sin(qJ(1));
t43 = cos(qJ(1));
t45 = pkin(4) + pkin(6);
t46 = pkin(4) - pkin(6);
t44 = cos(t46) / 0.2e1 + cos(t45) / 0.2e1;
t27 = t42 * t39 - t43 * t44;
t35 = -sin(t46) / 0.2e1;
t36 = sin(t45);
t34 = -t36 / 0.2e1 + t35;
t25 = atan2(-t27, t34);
t22 = sin(t25);
t23 = cos(t25);
t21 = -t22 * t27 + t23 * t34;
t29 = t43 * t39 + t42 * t44;
t47 = t29 ^ 2 / t21 ^ 2;
t41 = cos(pkin(6));
t40 = sin(pkin(4));
t38 = 0.1e1 / t42 ^ 2;
t33 = t36 / 0.2e1 + t35;
t32 = 0.1e1 / t34;
t30 = -t42 * t33 + t43 * t41;
t24 = 0.1e1 / (0.1e1 + t27 ^ 2 / t34 ^ 2);
t1 = [-t29 * t32 * t24, 0, 0, 0; (-t27 / t21 - (-t22 + (t23 * t27 * t32 + t22) * t24) * t47) / (0.1e1 + t47) 0, 0, 0; (-t41 + (-t33 / t42 - t30 * t38) * t43) / t40 / (0.1e1 + t30 ^ 2 * t38 / t40 ^ 2) 0, 0, 0;];
Ja_rot  = t1;
