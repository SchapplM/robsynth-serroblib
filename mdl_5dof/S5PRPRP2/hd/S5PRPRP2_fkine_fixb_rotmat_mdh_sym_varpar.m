% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5PRPRP2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:23
% EndTime: 2019-12-05 15:30:23
% DurationCPUTime: 0.13s
% Computational Cost: add. (118->37), mult. (83->38), div. (0->0), fcn. (125->8), ass. (0->32)
t17 = pkin(7) + qJ(2);
t12 = sin(t17);
t23 = sin(qJ(4));
t36 = t12 * t23;
t18 = sin(pkin(8));
t35 = t18 * t23;
t20 = cos(pkin(8));
t34 = t20 * t23;
t24 = cos(qJ(4));
t33 = t20 * t24;
t19 = sin(pkin(7));
t32 = t19 * pkin(1) + 0;
t21 = cos(pkin(7));
t31 = t21 * pkin(1) + 0;
t30 = qJ(1) + 0;
t29 = t12 * pkin(2) + t32;
t14 = pkin(5) + t30;
t13 = cos(t17);
t28 = t13 * pkin(2) + t12 * qJ(3) + t31;
t27 = pkin(3) * t20 + pkin(6) * t18;
t11 = t24 * pkin(4) + pkin(3);
t22 = -qJ(5) - pkin(6);
t26 = t11 * t20 - t18 * t22;
t25 = -t13 * qJ(3) + t29;
t10 = t18 * t24;
t6 = t13 * t18;
t5 = t12 * t18;
t4 = t13 * t33 + t36;
t3 = t12 * t24 - t13 * t34;
t2 = t12 * t33 - t13 * t23;
t1 = -t12 * t34 - t13 * t24;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t21, -t19, 0, 0; t19, t21, 0, 0; 0, 0, 1, t30; 0, 0, 0, 1; t13, -t12, 0, t31; t12, t13, 0, t32; 0, 0, 1, t14; 0, 0, 0, 1; t13 * t20, -t6, t12, t28; t12 * t20, -t5, -t13, t25; t18, t20, 0, t14; 0, 0, 0, 1; t4, t3, t6, t13 * t27 + t28; t2, t1, t5, t12 * t27 + t25; t10, -t35, -t20, t18 * pkin(3) - t20 * pkin(6) + t14; 0, 0, 0, 1; t4, t3, t6, pkin(4) * t36 + t13 * t26 + t28; t2, t1, t5, (-pkin(4) * t23 - qJ(3)) * t13 + t26 * t12 + t29; t10, -t35, -t20, t18 * t11 + t20 * t22 + t14; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
