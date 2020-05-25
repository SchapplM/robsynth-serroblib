% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPRRP6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:03
% EndTime: 2019-12-31 18:42:03
% DurationCPUTime: 0.12s
% Computational Cost: add. (118->37), mult. (83->38), div. (0->0), fcn. (125->8), ass. (0->32)
t17 = qJ(1) + pkin(8);
t12 = sin(t17);
t19 = sin(qJ(4));
t36 = t12 * t19;
t23 = cos(qJ(3));
t35 = t19 * t23;
t20 = sin(qJ(3));
t34 = t20 * t19;
t22 = cos(qJ(4));
t33 = t22 * t23;
t32 = pkin(5) + 0;
t21 = sin(qJ(1));
t31 = t21 * pkin(1) + 0;
t24 = cos(qJ(1));
t30 = t24 * pkin(1) + 0;
t29 = t12 * pkin(2) + t31;
t14 = qJ(2) + t32;
t13 = cos(t17);
t28 = t13 * pkin(2) + t12 * pkin(6) + t30;
t27 = pkin(3) * t23 + pkin(7) * t20;
t11 = t22 * pkin(4) + pkin(3);
t18 = -qJ(5) - pkin(7);
t26 = t11 * t23 - t18 * t20;
t25 = -t13 * pkin(6) + t29;
t10 = t20 * t22;
t6 = t13 * t20;
t5 = t12 * t20;
t4 = t13 * t33 + t36;
t3 = t12 * t22 - t13 * t35;
t2 = t12 * t33 - t13 * t19;
t1 = -t12 * t35 - t13 * t22;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t24, -t21, 0, 0; t21, t24, 0, 0; 0, 0, 1, t32; 0, 0, 0, 1; t13, -t12, 0, t30; t12, t13, 0, t31; 0, 0, 1, t14; 0, 0, 0, 1; t13 * t23, -t6, t12, t28; t12 * t23, -t5, -t13, t25; t20, t23, 0, t14; 0, 0, 0, 1; t4, t3, t6, t27 * t13 + t28; t2, t1, t5, t27 * t12 + t25; t10, -t34, -t23, t20 * pkin(3) - t23 * pkin(7) + t14; 0, 0, 0, 1; t4, t3, t6, pkin(4) * t36 + t26 * t13 + t28; t2, t1, t5, (-pkin(4) * t19 - pkin(6)) * t13 + t26 * t12 + t29; t10, -t34, -t23, t20 * t11 + t23 * t18 + t14; 0, 0, 0, 1;];
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
