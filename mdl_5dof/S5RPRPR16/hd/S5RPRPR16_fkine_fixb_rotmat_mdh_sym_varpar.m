% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPRPR16_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:37
% EndTime: 2019-12-31 18:38:37
% DurationCPUTime: 0.09s
% Computational Cost: add. (62->39), mult. (73->35), div. (0->0), fcn. (110->6), ass. (0->22)
t13 = sin(qJ(3));
t14 = sin(qJ(1));
t3 = t14 * t13;
t16 = cos(qJ(3));
t29 = t14 * t16;
t17 = cos(qJ(1));
t28 = t17 * t13;
t27 = t17 * t16;
t26 = qJ(4) * t16;
t11 = pkin(5) + 0;
t25 = t14 * pkin(1) + 0;
t24 = pkin(2) + t11;
t23 = t17 * pkin(1) + t14 * qJ(2) + 0;
t6 = t14 * pkin(6);
t22 = t17 * t26 + t25 + t6;
t21 = t17 * pkin(6) + t23;
t20 = t16 * pkin(3) + t13 * qJ(4) + t24;
t19 = pkin(3) * t3 + t21;
t18 = -t17 * qJ(2) + t25;
t15 = cos(qJ(5));
t12 = sin(qJ(5));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t17, -t14, 0, 0; t14, t17, 0, 0; 0, 0, 1, t11; 0, 0, 0, 1; 0, -t17, t14, t23; 0, -t14, -t17, t18; 1, 0, 0, t11; 0, 0, 0, 1; t3, t29, t17, t21; -t28, -t27, t14, t18 + t6; t16, -t13, 0, t24; 0, 0, 0, 1; t17, -t3, -t29, -t14 * t26 + t19; t14, t28, t27, (-pkin(3) * t13 - qJ(2)) * t17 + t22; 0, -t16, t13, t20; 0, 0, 0, 1; -t12 * t29 + t17 * t15, -t17 * t12 - t15 * t29, t3, t17 * pkin(4) + (pkin(7) * t13 - t26) * t14 + t19; t12 * t27 + t14 * t15, -t14 * t12 + t15 * t27, -t28, t14 * pkin(4) + (-qJ(2) + (-pkin(3) - pkin(7)) * t13) * t17 + t22; t13 * t12, t13 * t15, t16, t16 * pkin(7) + t20; 0, 0, 0, 1;];
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
