% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2018-11-23 15:53
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRPPR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:53:10
% EndTime: 2018-11-23 15:53:10
% DurationCPUTime: 0.10s
% Computational Cost: add. (160->44), mult. (108->40), div. (0->0), fcn. (153->8), ass. (0->31)
t22 = qJ(1) + pkin(9);
t15 = sin(t22);
t27 = cos(qJ(3));
t8 = t15 * t27;
t16 = cos(t22);
t10 = t16 * t27;
t23 = sin(qJ(6));
t24 = sin(qJ(3));
t43 = t23 * t24;
t26 = cos(qJ(6));
t42 = t24 * t26;
t41 = qJ(4) * t24;
t40 = pkin(6) + 0;
t25 = sin(qJ(1));
t39 = t25 * pkin(1) + 0;
t28 = cos(qJ(1));
t38 = t28 * pkin(1) + 0;
t17 = qJ(2) + t40;
t37 = t24 * pkin(3) + t17;
t36 = t16 * pkin(2) + t15 * pkin(7) + t38;
t35 = pkin(5) * t24 + pkin(8) * t27;
t34 = t15 * pkin(2) - t16 * pkin(7) + t39;
t33 = pkin(3) * t10 + t16 * t41 + t36;
t32 = -t27 * qJ(4) + t37;
t31 = pkin(3) * t8 + t15 * t41 + t34;
t30 = pkin(4) * t8 + t16 * qJ(5) + t31;
t29 = pkin(4) * t10 - t15 * qJ(5) + t33;
t18 = t24 * pkin(4);
t9 = t16 * t24;
t7 = t15 * t24;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t25, 0, 0; t25, t28, 0, 0; 0, 0, 1, t40; 0, 0, 0, 1; t16, -t15, 0, t38; t15, t16, 0, t39; 0, 0, 1, t17; 0, 0, 0, 1; t10, -t9, t15, t36; t8, -t7, -t16, t34; t24, t27, 0, t17; 0, 0, 0, 1; t10, t15, t9, t33; t8, -t16, t7, t31; t24, 0, -t27, t32; 0, 0, 0, 1; t9, -t10, -t15, t29; t7, -t8, t16, t30; -t27, -t24, 0, t18 + t32; 0, 0, 0, 1; -t15 * t23 + t16 * t42, -t15 * t26 - t16 * t43, t10, t35 * t16 + t29; t15 * t42 + t16 * t23, -t15 * t43 + t16 * t26, t8, t35 * t15 + t30; -t27 * t26, t27 * t23, t24, t24 * pkin(8) + t18 + (-pkin(5) - qJ(4)) * t27 + t37; 0, 0, 0, 1;];
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
