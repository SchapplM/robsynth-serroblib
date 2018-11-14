% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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

function Ja_rot = S4RPPP1_jacobia_rot_4_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobia_rot_4_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobia_rot_4_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:52
% EndTime: 2018-11-14 13:45:52
% DurationCPUTime: 0.05s
% Computational Cost: add. (199->19), mult. (196->29), div. (24->11), fcn. (226->12), ass. (0->24)
t43 = cos(pkin(6));
t44 = sin(qJ(1));
t45 = cos(qJ(1));
t47 = pkin(4) + pkin(6);
t48 = pkin(4) - pkin(6);
t46 = sin(t47) / 0.2e1 - sin(t48) / 0.2e1;
t29 = t44 * t43 + t45 * t46;
t37 = cos(t48) / 0.2e1;
t38 = cos(t47);
t36 = t37 - t38 / 0.2e1;
t27 = atan2(-t29, t36);
t24 = sin(t27);
t25 = cos(t27);
t23 = -t24 * t29 + t25 * t36;
t32 = t45 * t43 - t44 * t46;
t49 = t32 ^ 2 / t23 ^ 2;
t42 = sin(pkin(4));
t41 = sin(pkin(6));
t40 = 0.1e1 / t44 ^ 2;
t35 = t37 + t38 / 0.2e1;
t34 = 0.1e1 / t36;
t31 = -t44 * t35 - t45 * t41;
t26 = 0.1e1 / (0.1e1 + t29 ^ 2 / t36 ^ 2);
t1 = [-t32 * t34 * t26, 0, 0, 0; (-t29 / t23 - (-t24 + (t25 * t29 * t34 + t24) * t26) * t49) / (0.1e1 + t49) 0, 0, 0; (t41 + (-t35 / t44 - t31 * t40) * t45) / t42 / (0.1e1 + t31 ^ 2 * t40 / t42 ^ 2) 0, 0, 0;];
Ja_rot  = t1;
