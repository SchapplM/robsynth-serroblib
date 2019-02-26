% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
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

function Ja_rot = S4RPPP1_jacobia_rot_2_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobia_rot_2_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobia_rot_2_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:47
% EndTime: 2018-11-14 13:45:47
% DurationCPUTime: 0.05s
% Computational Cost: add. (52->17), mult. (98->36), div. (20->9), fcn. (140->13), ass. (0->23)
t31 = cos(pkin(4));
t29 = sin(pkin(4));
t33 = cos(qJ(1));
t34 = t33 * t29;
t23 = atan2(t34, t31);
t20 = sin(t23);
t21 = cos(t23);
t14 = t20 * t34 + t21 * t31;
t32 = sin(qJ(1));
t36 = 0.1e1 / t14 ^ 2 * t32 ^ 2;
t26 = t29 ^ 2;
t22 = 0.1e1 / (0.1e1 + t33 ^ 2 * t26 / t31 ^ 2);
t35 = t22 / t31;
t30 = cos(pkin(6));
t28 = sin(pkin(6));
t25 = pkin(4) - pkin(6);
t24 = pkin(4) + pkin(6);
t19 = cos(t25) / 0.2e1 + cos(t24) / 0.2e1;
t18 = sin(t24) / 0.2e1 - sin(t25) / 0.2e1;
t17 = -t32 * t18 + t33 * t30;
t16 = t32 * t19 + t33 * t28;
t15 = 0.1e1 / t17 ^ 2;
t1 = [-t32 * t29 * t35, 0, 0, 0; (0.1e1 / t14 * t34 - (-t21 * t26 * t33 * t35 + (t22 - 0.1e1) * t29 * t20) * t29 * t36) / (t26 * t36 + 0.1e1) 0, 0, 0; ((t33 * t19 - t32 * t28) / t17 - (-t33 * t18 - t32 * t30) * t16 * t15) / (t16 ^ 2 * t15 + 0.1e1) 0, 0, 0;];
Ja_rot  = t1;
