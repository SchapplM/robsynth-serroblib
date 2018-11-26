% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2018-11-23 16:24
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRRP2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:24:21
% EndTime: 2018-11-23 16:24:21
% DurationCPUTime: 0.15s
% Computational Cost: add. (199->55), mult. (127->60), div. (0->0), fcn. (182->10), ass. (0->42)
t31 = -pkin(9) - pkin(8);
t25 = sin(qJ(4));
t47 = t25 * pkin(4);
t28 = cos(qJ(4));
t13 = t28 * pkin(4) + pkin(3);
t23 = qJ(1) + pkin(10);
t14 = sin(t23);
t46 = t14 * t25;
t24 = qJ(4) + qJ(5);
t17 = sin(t24);
t29 = cos(qJ(3));
t45 = t17 * t29;
t18 = cos(t24);
t44 = t18 * t29;
t43 = t25 * t29;
t26 = sin(qJ(3));
t42 = t26 * t17;
t41 = t28 * t29;
t40 = pkin(6) + 0;
t27 = sin(qJ(1));
t39 = t27 * pkin(1) + 0;
t30 = cos(qJ(1));
t38 = t30 * pkin(1) + 0;
t16 = qJ(2) + t40;
t37 = t14 * pkin(2) + t39;
t15 = cos(t23);
t36 = t15 * pkin(2) + t14 * pkin(7) + t38;
t35 = pkin(3) * t29 + pkin(8) * t26;
t22 = -qJ(6) + t31;
t5 = pkin(5) * t18 + t13;
t34 = -t22 * t26 + t29 * t5;
t33 = t13 * t29 - t26 * t31;
t32 = -t15 * pkin(7) + t37;
t9 = t26 * t18;
t8 = t15 * t26;
t7 = t14 * t26;
t6 = pkin(5) * t17 + t47;
t4 = t14 * t17 + t15 * t44;
t3 = t14 * t18 - t15 * t45;
t2 = t14 * t44 - t15 * t17;
t1 = -t14 * t45 - t15 * t18;
t10 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t30, -t27, 0, 0; t27, t30, 0, 0; 0, 0, 1, t40; 0, 0, 0, 1; t15, -t14, 0, t38; t14, t15, 0, t39; 0, 0, 1, t16; 0, 0, 0, 1; t15 * t29, -t8, t14, t36; t14 * t29, -t7, -t15, t32; t26, t29, 0, t16; 0, 0, 0, 1; t15 * t41 + t46, t14 * t28 - t15 * t43, t8, t35 * t15 + t36; t14 * t41 - t15 * t25, -t14 * t43 - t15 * t28, t7, t35 * t14 + t32; t26 * t28, -t26 * t25, -t29, t26 * pkin(3) - t29 * pkin(8) + t16; 0, 0, 0, 1; t4, t3, t8, pkin(4) * t46 + t33 * t15 + t36; t2, t1, t7 (-pkin(7) - t47) * t15 + t33 * t14 + t37; t9, -t42, -t29, t26 * t13 + t29 * t31 + t16; 0, 0, 0, 1; t4, t3, t8, t14 * t6 + t34 * t15 + t36; t2, t1, t7 (-pkin(7) - t6) * t15 + t34 * t14 + t37; t9, -t42, -t29, t29 * t22 + t26 * t5 + t16; 0, 0, 0, 1;];
T_ges = t10;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
