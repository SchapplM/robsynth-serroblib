% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2018-11-23 17:00
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPRPR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:00:23
% EndTime: 2018-11-23 17:00:23
% DurationCPUTime: 0.14s
% Computational Cost: add. (186->55), mult. (100->50), div. (0->0), fcn. (149->10), ass. (0->35)
t30 = cos(qJ(1));
t22 = qJ(2) + pkin(10);
t18 = qJ(4) + t22;
t13 = sin(t18);
t38 = qJ(5) * t13;
t14 = cos(t18);
t9 = t30 * t14;
t45 = pkin(4) * t9 + t30 * t38;
t29 = cos(qJ(2));
t15 = t29 * pkin(2) + pkin(1);
t27 = sin(qJ(1));
t44 = t27 * t13;
t8 = t27 * t14;
t25 = sin(qJ(6));
t43 = t27 * t25;
t28 = cos(qJ(6));
t42 = t27 * t28;
t41 = t30 * t13;
t40 = t30 * t25;
t39 = t30 * t28;
t24 = -qJ(3) - pkin(7);
t23 = pkin(6) + 0;
t17 = cos(t22);
t3 = pkin(3) * t17 + t15;
t37 = t30 * t3 + 0;
t21 = -pkin(8) + t24;
t36 = t30 * t21 + t27 * t3 + 0;
t26 = sin(qJ(2));
t35 = t26 * pkin(2) + t23;
t16 = sin(t22);
t34 = pkin(3) * t16 + t35;
t33 = pkin(4) * t8 + t27 * t38 + t36;
t32 = -t27 * t21 + t37;
t31 = t13 * pkin(4) - t14 * qJ(5) + t34;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t30, -t27, 0, 0; t27, t30, 0, 0; 0, 0, 1, t23; 0, 0, 0, 1; t30 * t29, -t30 * t26, t27, t30 * pkin(1) + t27 * pkin(7) + 0; t27 * t29, -t27 * t26, -t30, t27 * pkin(1) - t30 * pkin(7) + 0; t26, t29, 0, t23; 0, 0, 0, 1; t30 * t17, -t30 * t16, t27, t30 * t15 - t27 * t24 + 0; t27 * t17, -t27 * t16, -t30, t27 * t15 + t30 * t24 + 0; t16, t17, 0, t35; 0, 0, 0, 1; t9, -t41, t27, t32; t8, -t44, -t30, t36; t13, t14, 0, t34; 0, 0, 0, 1; t27, -t9, t41, t32 + t45; -t30, -t8, t44, t33; 0, -t13, -t14, t31; 0, 0, 0, 1; t13 * t40 + t42, t13 * t39 - t43, t9, pkin(9) * t9 + (pkin(5) - t21) * t27 + t37 + t45; t13 * t43 - t39, t13 * t42 + t40, t8, -t30 * pkin(5) + pkin(9) * t8 + t33; -t14 * t25, -t14 * t28, t13, t13 * pkin(9) + t31; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
