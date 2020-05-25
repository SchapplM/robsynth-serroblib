% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2018-11-23 15:54
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRPPR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:53:39
% EndTime: 2018-11-23 15:53:40
% DurationCPUTime: 0.13s
% Computational Cost: add. (201->57), mult. (199->61), div. (0->0), fcn. (278->10), ass. (0->40)
t25 = pkin(9) + qJ(3);
t22 = sin(t25);
t29 = cos(pkin(10));
t15 = t22 * t29;
t27 = sin(pkin(10));
t52 = t22 * t27;
t53 = pkin(4) * t15 + qJ(5) * t52;
t33 = sin(qJ(1));
t16 = t33 * t22;
t23 = cos(t25);
t51 = t33 * t23;
t50 = t33 * t27;
t49 = t33 * t29;
t35 = cos(qJ(1));
t17 = t35 * t22;
t48 = t35 * t23;
t47 = t35 * t27;
t46 = t35 * t29;
t45 = qJ(4) * t22;
t26 = pkin(6) + 0;
t28 = sin(pkin(9));
t44 = t28 * pkin(2) + t26;
t30 = cos(pkin(9));
t20 = t30 * pkin(2) + pkin(1);
t31 = -pkin(7) - qJ(2);
t43 = t33 * t20 + t35 * t31 + 0;
t42 = t22 * pkin(3) + t44;
t41 = t35 * t20 - t33 * t31 + 0;
t40 = pkin(3) * t51 + t33 * t45 + t43;
t39 = -t23 * qJ(4) + t42;
t38 = pkin(3) * t48 + t35 * t45 + t41;
t3 = t23 * t50 + t46;
t4 = t23 * t49 - t47;
t37 = t4 * pkin(4) + t3 * qJ(5) + t40;
t5 = t23 * t47 - t49;
t6 = t23 * t46 + t50;
t36 = t6 * pkin(4) + t5 * qJ(5) + t38;
t34 = cos(qJ(6));
t32 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t35, -t33, 0, 0; t33, t35, 0, 0; 0, 0, 1, t26; 0, 0, 0, 1; t35 * t30, -t35 * t28, t33, t35 * pkin(1) + t33 * qJ(2) + 0; t33 * t30, -t33 * t28, -t35, t33 * pkin(1) - t35 * qJ(2) + 0; t28, t30, 0, t26; 0, 0, 0, 1; t48, -t17, t33, t41; t51, -t16, -t35, t43; t22, t23, 0, t44; 0, 0, 0, 1; t6, -t5, t17, t38; t4, -t3, t16, t40; t15, -t52, -t23, t39; 0, 0, 0, 1; t6, t17, t5, t36; t4, t16, t3, t37; t15, -t23, t52, t39 + t53; 0, 0, 0, 1; t5 * t32 + t6 * t34, -t6 * t32 + t5 * t34, -t17, t6 * pkin(5) - pkin(8) * t17 + t36; t3 * t32 + t4 * t34, t3 * t34 - t4 * t32, -t16, t4 * pkin(5) - pkin(8) * t16 + t37; (t27 * t32 + t29 * t34) * t22 (t27 * t34 - t29 * t32) * t22, t23, pkin(5) * t15 + (pkin(8) - qJ(4)) * t23 + t42 + t53; 0, 0, 0, 1;];
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
