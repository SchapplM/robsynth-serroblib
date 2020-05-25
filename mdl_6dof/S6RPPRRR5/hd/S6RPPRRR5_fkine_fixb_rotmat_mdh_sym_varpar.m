% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:49
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPPRRR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:49:32
% EndTime: 2018-11-23 15:49:32
% DurationCPUTime: 0.10s
% Computational Cost: add. (100->42), mult. (78->37), div. (0->0), fcn. (119->8), ass. (0->32)
t16 = sin(qJ(1));
t13 = qJ(4) + qJ(5);
t5 = cos(t13);
t39 = t16 * t5;
t19 = cos(qJ(1));
t38 = t19 * t5;
t14 = sin(qJ(6));
t37 = t16 * t14;
t15 = sin(qJ(4));
t36 = t16 * t15;
t17 = cos(qJ(6));
t35 = t16 * t17;
t34 = t19 * t14;
t33 = t19 * t15;
t32 = t19 * t17;
t12 = pkin(6) + 0;
t31 = t16 * pkin(1) + 0;
t30 = pkin(2) + t12;
t6 = t16 * qJ(3);
t29 = t6 + t31;
t28 = t19 * pkin(1) + t16 * qJ(2) + 0;
t27 = pkin(3) + t30;
t26 = t19 * qJ(3) + t28;
t4 = sin(t13);
t25 = pkin(5) * t4 - pkin(9) * t5;
t18 = cos(qJ(4));
t24 = t18 * pkin(4) + t27;
t23 = -t19 * qJ(2) + t31;
t20 = -pkin(8) - pkin(7);
t22 = pkin(4) * t33 + t16 * t20 + t26;
t21 = pkin(4) * t36 + (-qJ(2) - t20) * t19 + t29;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t19, -t16, 0, 0; t16, t19, 0, 0; 0, 0, 1, t12; 0, 0, 0, 1; 0, -t19, t16, t28; 0, -t16, -t19, t23; 1, 0, 0, t12; 0, 0, 0, 1; 0, t16, t19, t26; 0, -t19, t16, t23 + t6; 1, 0, 0, t30; 0, 0, 0, 1; t33, t19 * t18, -t16, -t16 * pkin(7) + t26; t36, t16 * t18, t19 (pkin(7) - qJ(2)) * t19 + t29; t18, -t15, 0, t27; 0, 0, 0, 1; t19 * t4, t38, -t16, t22; t16 * t4, t39, t19, t21; t5, -t4, 0, t24; 0, 0, 0, 1; t4 * t32 - t37, -t4 * t34 - t35, -t38, t25 * t19 + t22; t4 * t35 + t34, -t4 * t37 + t32, -t39, t25 * t16 + t21; t5 * t17, -t5 * t14, t4, t5 * pkin(5) + t4 * pkin(9) + t24; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
