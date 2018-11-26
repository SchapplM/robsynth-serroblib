% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2018-11-23 16:05
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRPRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:04:35
% EndTime: 2018-11-23 16:04:35
% DurationCPUTime: 0.17s
% Computational Cost: add. (179->56), mult. (120->56), div. (0->0), fcn. (170->10), ass. (0->37)
t19 = qJ(1) + pkin(10);
t11 = sin(t19);
t22 = sin(qJ(3));
t39 = qJ(4) * t22;
t25 = cos(qJ(3));
t5 = t11 * t25;
t46 = pkin(3) * t5 + t11 * t39;
t45 = t11 * t22;
t12 = cos(t19);
t44 = t12 * t22;
t6 = t12 * t25;
t20 = qJ(5) + qJ(6);
t14 = sin(t20);
t43 = t14 * t22;
t15 = cos(t20);
t42 = t15 * t22;
t21 = sin(qJ(5));
t41 = t21 * t22;
t24 = cos(qJ(5));
t40 = t22 * t24;
t38 = pkin(6) + 0;
t23 = sin(qJ(1));
t37 = t23 * pkin(1) + 0;
t26 = cos(qJ(1));
t36 = t26 * pkin(1) + 0;
t35 = t11 * pkin(2) + t37;
t13 = qJ(2) + t38;
t34 = t12 * pkin(2) + t11 * pkin(7) + t36;
t33 = t22 * pkin(3) + t13;
t32 = t35 + t46;
t31 = pkin(3) * t6 + t12 * t39 + t34;
t27 = -pkin(9) - pkin(8);
t30 = pkin(5) * t41 - t25 * t27;
t29 = -t12 * pkin(7) + t35;
t28 = -t25 * qJ(4) + t33;
t10 = t24 * pkin(5) + pkin(4);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t26, -t23, 0, 0; t23, t26, 0, 0; 0, 0, 1, t38; 0, 0, 0, 1; t12, -t11, 0, t36; t11, t12, 0, t37; 0, 0, 1, t13; 0, 0, 0, 1; t6, -t44, t11, t34; t5, -t45, -t12, t29; t22, t25, 0, t13; 0, 0, 0, 1; t11, -t6, t44, t31; -t12, -t5, t45, t29 + t46; 0, -t22, -t25, t28; 0, 0, 0, 1; t11 * t24 + t12 * t41, -t11 * t21 + t12 * t40, t6, t11 * pkin(4) + pkin(8) * t6 + t31; t11 * t41 - t12 * t24, t11 * t40 + t12 * t21, t5, pkin(8) * t5 + (-pkin(4) - pkin(7)) * t12 + t32; -t25 * t21, -t25 * t24, t22, t22 * pkin(8) + t28; 0, 0, 0, 1; t11 * t15 + t12 * t43, -t11 * t14 + t12 * t42, t6, t11 * t10 + t30 * t12 + t31; t11 * t43 - t12 * t15, t11 * t42 + t12 * t14, t5 (-pkin(7) - t10) * t12 + t30 * t11 + t32; -t25 * t14, -t25 * t15, t22, -t22 * t27 + (-pkin(5) * t21 - qJ(4)) * t25 + t33; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
