% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2018-11-23 15:44
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPPRRP3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:44:40
% EndTime: 2018-11-23 15:44:40
% DurationCPUTime: 0.10s
% Computational Cost: add. (159->43), mult. (108->38), div. (0->0), fcn. (158->8), ass. (0->34)
t22 = qJ(1) + pkin(9);
t15 = sin(t22);
t24 = sin(qJ(4));
t44 = t15 * t24;
t27 = cos(qJ(4));
t43 = t15 * t27;
t16 = cos(t22);
t8 = t16 * t27;
t23 = sin(qJ(5));
t42 = t23 * t24;
t26 = cos(qJ(5));
t41 = t24 * t26;
t40 = t27 * t23;
t39 = pkin(6) + 0;
t25 = sin(qJ(1));
t38 = t25 * pkin(1) + 0;
t28 = cos(qJ(1));
t37 = t28 * pkin(1) + 0;
t17 = qJ(2) + t39;
t36 = t15 * pkin(2) + t38;
t35 = pkin(3) + t17;
t34 = t16 * pkin(2) + t15 * qJ(3) + t37;
t33 = t16 * pkin(7) + t34;
t32 = t27 * pkin(4) + t24 * pkin(8) + t35;
t31 = -t16 * qJ(3) + t36;
t30 = pkin(4) * t44 - pkin(8) * t43 + t33;
t10 = t15 * pkin(7);
t29 = t10 + pkin(8) * t8 + (-pkin(4) * t24 - qJ(3)) * t16 + t36;
t14 = t27 * t26;
t4 = t15 * t23 - t16 * t41;
t3 = t15 * t26 + t16 * t42;
t2 = t15 * t41 + t16 * t23;
t1 = t15 * t42 - t16 * t26;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t25, 0, 0; t25, t28, 0, 0; 0, 0, 1, t39; 0, 0, 0, 1; t16, -t15, 0, t37; t15, t16, 0, t38; 0, 0, 1, t17; 0, 0, 0, 1; 0, -t16, t15, t34; 0, -t15, -t16, t31; 1, 0, 0, t17; 0, 0, 0, 1; t44, t43, t16, t33; -t16 * t24, -t8, t15, t10 + t31; t27, -t24, 0, t35; 0, 0, 0, 1; t2, -t1, -t43, t30; t4, t3, t8, t29; t14, -t40, t24, t32; 0, 0, 0, 1; t2, -t43, t1, t2 * pkin(5) + t1 * qJ(6) + t30; t4, t8, -t3, t4 * pkin(5) - t3 * qJ(6) + t29; t14, t24, t40 (pkin(5) * t26 + qJ(6) * t23) * t27 + t32; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
