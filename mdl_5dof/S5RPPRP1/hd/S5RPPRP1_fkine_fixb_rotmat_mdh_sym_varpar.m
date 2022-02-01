% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% Tc_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)
% T_c_stack [(5+1)*3 x 4]
%   stacked matrices from Tc_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RPPRP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:05
% EndTime: 2022-01-23 09:12:05
% DurationCPUTime: 0.14s
% Computational Cost: add. (118->37), mult. (83->38), div. (0->0), fcn. (125->8), ass. (0->32)
t17 = qJ(1) + pkin(7);
t12 = sin(t17);
t21 = sin(qJ(4));
t36 = t12 * t21;
t18 = sin(pkin(8));
t35 = t18 * t21;
t19 = cos(pkin(8));
t34 = t19 * t21;
t23 = cos(qJ(4));
t33 = t19 * t23;
t32 = pkin(5) + 0;
t22 = sin(qJ(1));
t31 = t22 * pkin(1) + 0;
t24 = cos(qJ(1));
t30 = t24 * pkin(1) + 0;
t29 = t12 * pkin(2) + t31;
t14 = qJ(2) + t32;
t13 = cos(t17);
t28 = t13 * pkin(2) + t12 * qJ(3) + t30;
t27 = pkin(3) * t19 + pkin(6) * t18;
t11 = t23 * pkin(4) + pkin(3);
t20 = -qJ(5) - pkin(6);
t26 = t11 * t19 - t18 * t20;
t25 = -t13 * qJ(3) + t29;
t10 = t18 * t23;
t6 = t13 * t18;
t5 = t12 * t18;
t4 = t13 * t33 + t36;
t3 = t12 * t23 - t13 * t34;
t2 = t12 * t33 - t13 * t21;
t1 = -t12 * t34 - t13 * t23;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t24, -t22, 0, 0; t22, t24, 0, 0; 0, 0, 1, t32; t13, -t12, 0, t30; t12, t13, 0, t31; 0, 0, 1, t14; t13 * t19, -t6, t12, t28; t12 * t19, -t5, -t13, t25; t18, t19, 0, t14; t4, t3, t6, t13 * t27 + t28; t2, t1, t5, t12 * t27 + t25; t10, -t35, -t19, t18 * pkin(3) - t19 * pkin(6) + t14; t4, t3, t6, pkin(4) * t36 + t13 * t26 + t28; t2, t1, t5, (-pkin(4) * t21 - qJ(3)) * t13 + t26 * t12 + t29; t10, -t35, -t19, t18 * t11 + t19 * t20 + t14;];
Tc_stack = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
