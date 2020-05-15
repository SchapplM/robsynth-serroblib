% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPRRP12_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:09
% EndTime: 2019-12-31 18:56:09
% DurationCPUTime: 0.12s
% Computational Cost: add. (69->40), mult. (85->35), div. (0->0), fcn. (127->6), ass. (0->33)
t15 = sin(qJ(4));
t17 = sin(qJ(1));
t36 = t17 * t15;
t16 = sin(qJ(3));
t35 = t17 * t16;
t18 = cos(qJ(4));
t34 = t17 * t18;
t19 = cos(qJ(3));
t33 = t17 * t19;
t32 = t19 * t15;
t20 = cos(qJ(1));
t31 = t20 * t15;
t30 = t20 * t16;
t29 = t20 * t18;
t13 = pkin(5) + 0;
t28 = t17 * pkin(1) + 0;
t27 = pkin(2) + t13;
t9 = t17 * pkin(6);
t26 = t9 + t28;
t25 = t20 * pkin(1) + t17 * qJ(2) + 0;
t24 = t20 * pkin(6) + t25;
t23 = pkin(3) * t16 - pkin(7) * t19;
t14 = -qJ(5) - pkin(7);
t7 = t18 * pkin(4) + pkin(3);
t22 = t14 * t19 + t16 * t7;
t21 = -t20 * qJ(2) + t28;
t6 = t20 * t19;
t5 = t19 * t18;
t4 = -t16 * t29 + t36;
t3 = t15 * t30 + t34;
t2 = t16 * t34 + t31;
t1 = -t15 * t35 + t29;
t8 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t20, -t17, 0, 0; t17, t20, 0, 0; 0, 0, 1, t13; 0, 0, 0, 1; 0, -t20, t17, t25; 0, -t17, -t20, t21; 1, 0, 0, t13; 0, 0, 0, 1; t35, t33, t20, t24; -t30, -t6, t17, t21 + t9; t19, -t16, 0, t27; 0, 0, 0, 1; t2, t1, -t33, t23 * t17 + t24; t4, t3, t6, (-qJ(2) - t23) * t20 + t26; t5, -t32, t16, t19 * pkin(3) + t16 * pkin(7) + t27; 0, 0, 0, 1; t2, t1, -t33, pkin(4) * t31 + t22 * t17 + t24; t4, t3, t6, pkin(4) * t36 + (-qJ(2) - t22) * t20 + t26; t5, -t32, t16, -t16 * t14 + t19 * t7 + t27; 0, 0, 0, 1;];
T_ges = t8;
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
