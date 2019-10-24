% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-10-24 10:24
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRPRP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:24:45
% EndTime: 2019-10-24 10:24:46
% DurationCPUTime: 0.14s
% Computational Cost: add. (116->47), mult. (132->52), div. (0->0), fcn. (187->8), ass. (0->32)
t22 = sin(pkin(8));
t23 = sin(pkin(7));
t40 = t23 * t22;
t28 = cos(qJ(2));
t39 = t23 * t28;
t25 = cos(pkin(7));
t38 = t25 * t28;
t26 = -pkin(6) - qJ(3);
t27 = sin(qJ(2));
t37 = t26 * t27;
t21 = pkin(8) + qJ(4);
t15 = sin(t21);
t36 = t27 * t15;
t35 = t23 * pkin(1) + 0;
t20 = qJ(1) + 0;
t34 = t25 * pkin(1) + t23 * pkin(5) + 0;
t24 = cos(pkin(8));
t13 = t24 * pkin(3) + pkin(2);
t33 = t27 * t13 + t28 * t26 + t20;
t32 = pkin(2) * t28 + qJ(3) * t27;
t31 = -t25 * pkin(5) + t35;
t30 = pkin(3) * t40 + t13 * t38 - t25 * t37 + t34;
t29 = -t23 * t37 + t13 * t39 + (-pkin(3) * t22 - pkin(5)) * t25 + t35;
t16 = cos(t21);
t12 = t25 * t27;
t11 = t23 * t27;
t9 = t27 * t16;
t4 = t23 * t15 + t16 * t38;
t3 = t15 * t38 - t23 * t16;
t2 = -t25 * t15 + t16 * t39;
t1 = t15 * t39 + t25 * t16;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t25, -t23, 0, 0; t23, t25, 0, 0; 0, 0, 1, t20; 0, 0, 0, 1; t38, -t12, t23, t34; t39, -t11, -t25, t31; t27, t28, 0, t20; 0, 0, 0, 1; t24 * t38 + t40, -t22 * t38 + t23 * t24, t12, t32 * t25 + t34; -t25 * t22 + t24 * t39, -t22 * t39 - t25 * t24, t11, t32 * t23 + t31; t27 * t24, -t27 * t22, -t28, t27 * pkin(2) - t28 * qJ(3) + t20; 0, 0, 0, 1; t4, -t3, t12, t30; t2, -t1, t11, t29; t9, -t36, -t28, t33; 0, 0, 0, 1; t4, t12, t3, t4 * pkin(4) + t3 * qJ(5) + t30; t2, t11, t1, t2 * pkin(4) + t1 * qJ(5) + t29; t9, -t28, t36, (pkin(4) * t16 + qJ(5) * t15) * t27 + t33; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
