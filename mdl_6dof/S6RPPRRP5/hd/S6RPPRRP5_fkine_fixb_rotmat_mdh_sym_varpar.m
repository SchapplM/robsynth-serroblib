% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2018-11-23 15:46
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPPRRP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:45:50
% EndTime: 2018-11-23 15:45:50
% DurationCPUTime: 0.10s
% Computational Cost: add. (88->46), mult. (97->38), div. (0->0), fcn. (143->6), ass. (0->36)
t15 = sin(qJ(5));
t39 = pkin(5) * t15;
t16 = sin(qJ(4));
t17 = sin(qJ(1));
t38 = t17 * t16;
t18 = cos(qJ(5));
t37 = t17 * t18;
t19 = cos(qJ(4));
t36 = t17 * t19;
t35 = t19 * t15;
t20 = cos(qJ(1));
t34 = t20 * t16;
t33 = t20 * t18;
t32 = t20 * t19;
t13 = pkin(6) + 0;
t31 = t17 * pkin(1) + 0;
t30 = pkin(2) + t13;
t29 = t20 * pkin(1) + t17 * qJ(2) + 0;
t28 = pkin(3) + t30;
t27 = t20 * qJ(3) + t29;
t26 = pkin(4) * t16 - pkin(8) * t19;
t14 = -qJ(6) - pkin(8);
t6 = t18 * pkin(5) + pkin(4);
t25 = t14 * t19 + t16 * t6;
t24 = -t20 * qJ(2) + t31;
t7 = t17 * qJ(3);
t23 = t24 + t7;
t22 = -t17 * pkin(7) + t27;
t11 = t20 * pkin(7);
t21 = t11 + t23;
t5 = t19 * t18;
t4 = -t17 * t15 + t16 * t33;
t3 = -t15 * t34 - t37;
t2 = t20 * t15 + t16 * t37;
t1 = -t15 * t38 + t33;
t8 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t20, -t17, 0, 0; t17, t20, 0, 0; 0, 0, 1, t13; 0, 0, 0, 1; 0, -t20, t17, t29; 0, -t17, -t20, t24; 1, 0, 0, t13; 0, 0, 0, 1; 0, t17, t20, t27; 0, -t20, t17, t23; 1, 0, 0, t30; 0, 0, 0, 1; t34, t32, -t17, t22; t38, t36, t20, t21; t19, -t16, 0, t28; 0, 0, 0, 1; t4, t3, -t32, t26 * t20 + t22; t2, t1, -t36, t26 * t17 + t21; t5, -t35, t16, t19 * pkin(4) + t16 * pkin(8) + t28; 0, 0, 0, 1; t4, t3, -t32, t25 * t20 + (-pkin(7) - t39) * t17 + t27; t2, t1, -t36, t11 + t7 + (-qJ(2) + t39) * t20 + t25 * t17 + t31; t5, -t35, t16, -t16 * t14 + t19 * t6 + t28; 0, 0, 0, 1;];
T_ges = t8;
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
