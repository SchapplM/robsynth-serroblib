% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% Datum: 2018-11-23 15:38
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPPPRR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:38:04
% EndTime: 2018-11-23 15:38:04
% DurationCPUTime: 0.13s
% Computational Cost: add. (144->49), mult. (164->46), div. (0->0), fcn. (247->10), ass. (0->30)
t19 = pkin(10) + qJ(5);
t11 = sin(t19);
t34 = sin(pkin(9));
t35 = cos(pkin(9));
t38 = sin(qJ(1));
t39 = cos(qJ(1));
t5 = -t34 * t38 - t35 * t39;
t41 = t5 * t11;
t6 = t34 * t39 - t35 * t38;
t40 = t6 * t11;
t12 = cos(t19);
t24 = sin(qJ(6));
t37 = t12 * t24;
t25 = cos(qJ(6));
t36 = t12 * t25;
t20 = pkin(6) + 0;
t13 = -qJ(3) + t20;
t33 = t39 * pkin(1) + t38 * qJ(2) + 0;
t32 = t39 * pkin(2) + t33;
t31 = -pkin(5) * t12 - pkin(8) * t11;
t21 = sin(pkin(10));
t30 = -t21 * pkin(4) + t13;
t22 = cos(pkin(10));
t10 = t22 * pkin(4) + pkin(3);
t23 = -pkin(7) - qJ(4);
t29 = -t5 * t10 - t6 * t23 + t32;
t28 = t38 * pkin(1) - qJ(2) * t39 + 0;
t27 = t38 * pkin(2) + t28;
t26 = -t6 * t10 + t5 * t23 + t27;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t39, -t38, 0, 0; t38, t39, 0, 0; 0, 0, 1, t20; 0, 0, 0, 1; t39, 0, t38, t33; t38, 0, -t39, t28; 0, 1, 0, t20; 0, 0, 0, 1; -t5, -t6, 0, t32; -t6, t5, 0, t27; 0, 0, -1, t13; 0, 0, 0, 1; -t5 * t22, t5 * t21, t6, -t5 * pkin(3) + t6 * qJ(4) + t32; -t6 * t22, t6 * t21, -t5, -t6 * pkin(3) - t5 * qJ(4) + t27; -t21, -t22, 0, t13; 0, 0, 0, 1; -t5 * t12, t41, t6, t29; -t6 * t12, t40, -t5, t26; -t11, -t12, 0, t30; 0, 0, 0, 1; t6 * t24 - t36 * t5, t6 * t25 + t37 * t5, -t41, t31 * t5 + t29; -t5 * t24 - t36 * t6, -t5 * t25 + t37 * t6, -t40, t31 * t6 + t26; -t11 * t25, t11 * t24, t12, -t11 * pkin(5) + t12 * pkin(8) + t30; 0, 0, 0, 1;];
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
