% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRPPR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:46
% EndTime: 2019-12-31 19:28:46
% DurationCPUTime: 0.10s
% Computational Cost: add. (108->39), mult. (84->38), div. (0->0), fcn. (126->8), ass. (0->27)
t23 = cos(qJ(1));
t15 = qJ(2) + pkin(8);
t12 = sin(t15);
t32 = qJ(4) * t12;
t13 = cos(t15);
t8 = t23 * t13;
t35 = pkin(3) * t8 + t23 * t32;
t20 = sin(qJ(1));
t34 = t20 * t12;
t7 = t20 * t13;
t33 = t23 * t12;
t16 = pkin(5) + 0;
t22 = cos(qJ(2));
t11 = t22 * pkin(2) + pkin(1);
t31 = t23 * t11 + 0;
t19 = sin(qJ(2));
t30 = t19 * pkin(2) + t16;
t17 = -qJ(3) - pkin(6);
t29 = t20 * t11 + t23 * t17 + 0;
t28 = pkin(3) * t7 + t20 * t32 + t29;
t18 = sin(qJ(5));
t21 = cos(qJ(5));
t27 = t12 * t21 - t13 * t18;
t26 = t12 * t18 + t13 * t21;
t25 = -t20 * t17 + t31;
t24 = t12 * pkin(3) - t13 * qJ(4) + t30;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t23, -t20, 0, 0; t20, t23, 0, 0; 0, 0, 1, t16; 0, 0, 0, 1; t23 * t22, -t23 * t19, t20, t23 * pkin(1) + t20 * pkin(6) + 0; t20 * t22, -t20 * t19, -t23, t20 * pkin(1) - t23 * pkin(6) + 0; t19, t22, 0, t16; 0, 0, 0, 1; t8, -t33, t20, t25; t7, -t34, -t23, t29; t12, t13, 0, t30; 0, 0, 0, 1; t8, t20, t33, t25 + t35; t7, -t23, t34, t28; t12, 0, -t13, t24; 0, 0, 0, 1; t26 * t23, t27 * t23, -t20, pkin(4) * t8 + (-pkin(7) - t17) * t20 + t31 + t35; t26 * t20, t27 * t20, t23, pkin(4) * t7 + t23 * pkin(7) + t28; t27, -t26, 0, t12 * pkin(4) + t24; 0, 0, 0, 1;];
T_ges = t1;
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
