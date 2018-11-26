% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2018-11-23 16:55
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPPRR10_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:55:08
% EndTime: 2018-11-23 16:55:08
% DurationCPUTime: 0.16s
% Computational Cost: add. (151->70), mult. (162->71), div. (0->0), fcn. (221->10), ass. (0->36)
t27 = sin(qJ(1));
t28 = cos(qJ(2));
t10 = t27 * t28;
t26 = sin(qJ(2));
t40 = qJ(3) * t26;
t44 = pkin(2) * t10 + t27 * t40;
t23 = sin(pkin(10));
t43 = t23 * pkin(4);
t24 = cos(pkin(10));
t9 = t24 * pkin(4) + pkin(3);
t42 = t27 * t26;
t29 = cos(qJ(1));
t41 = t29 * t26;
t11 = t29 * t28;
t25 = -pkin(8) - qJ(4);
t39 = qJ(4) * t28;
t22 = pkin(6) + 0;
t21 = pkin(10) + qJ(5);
t38 = t27 * pkin(1) + 0;
t37 = t26 * pkin(2) + t22;
t36 = t29 * pkin(1) + t27 * pkin(7) + 0;
t35 = t38 + t44;
t12 = sin(t21);
t2 = pkin(5) * t12 + t43;
t20 = -pkin(9) + t25;
t34 = t2 * t26 - t20 * t28;
t33 = -t29 * pkin(7) + t38;
t32 = pkin(2) * t11 + t29 * t40 + t36;
t31 = -t25 * t28 + t26 * t43;
t30 = -t28 * qJ(3) + t37;
t14 = qJ(6) + t21;
t13 = cos(t21);
t8 = cos(t14);
t7 = sin(t14);
t1 = pkin(5) * t13 + t9;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t27, 0, 0; t27, t29, 0, 0; 0, 0, 1, t22; 0, 0, 0, 1; t11, -t41, t27, t36; t10, -t42, -t29, t33; t26, t28, 0, t22; 0, 0, 0, 1; t27, -t11, t41, t32; -t29, -t10, t42, t33 + t44; 0, -t26, -t28, t30; 0, 0, 0, 1; t23 * t41 + t27 * t24, -t27 * t23 + t24 * t41, t11, t27 * pkin(3) + t29 * t39 + t32; t23 * t42 - t29 * t24, t29 * t23 + t24 * t42, t10, t27 * t39 + (-pkin(3) - pkin(7)) * t29 + t35; -t28 * t23, -t28 * t24, t26, t26 * qJ(4) + t30; 0, 0, 0, 1; t12 * t41 + t27 * t13, -t27 * t12 + t13 * t41, t11, t27 * t9 + t31 * t29 + t32; t12 * t42 - t29 * t13, t29 * t12 + t13 * t42, t10 (-pkin(7) - t9) * t29 + t31 * t27 + t35; -t28 * t12, -t28 * t13, t26, -t26 * t25 + (-qJ(3) - t43) * t28 + t37; 0, 0, 0, 1; t27 * t8 + t7 * t41, -t27 * t7 + t8 * t41, t11, t27 * t1 + t34 * t29 + t32; -t29 * t8 + t7 * t42, t29 * t7 + t8 * t42, t10 (-pkin(7) - t1) * t29 + t34 * t27 + t35; -t28 * t7, -t28 * t8, t26, -t26 * t20 + (-qJ(3) - t2) * t28 + t37; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
