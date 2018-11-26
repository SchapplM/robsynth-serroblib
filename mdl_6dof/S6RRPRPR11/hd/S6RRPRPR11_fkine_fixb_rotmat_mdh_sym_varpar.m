% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2018-11-23 17:08
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPRPR11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:07:43
% EndTime: 2018-11-23 17:07:43
% DurationCPUTime: 0.16s
% Computational Cost: add. (151->70), mult. (162->70), div. (0->0), fcn. (221->10), ass. (0->37)
t26 = sin(qJ(1));
t25 = sin(qJ(2));
t39 = qJ(3) * t25;
t28 = cos(qJ(2));
t9 = t26 * t28;
t45 = pkin(2) * t9 + t26 * t39;
t24 = sin(qJ(4));
t44 = t24 * pkin(4);
t27 = cos(qJ(4));
t11 = t27 * pkin(4) + pkin(3);
t43 = t26 * t25;
t42 = t26 * t27;
t29 = cos(qJ(1));
t41 = t29 * t25;
t40 = t29 * t27;
t10 = t29 * t28;
t23 = -qJ(5) - pkin(8);
t22 = pkin(6) + 0;
t21 = qJ(4) + pkin(10);
t38 = t26 * pkin(1) + 0;
t37 = t25 * pkin(2) + t22;
t36 = t29 * pkin(1) + t26 * pkin(7) + 0;
t35 = t38 + t45;
t12 = sin(t21);
t2 = pkin(5) * t12 + t44;
t20 = -pkin(9) + t23;
t34 = t2 * t25 - t20 * t28;
t33 = -t29 * pkin(7) + t38;
t32 = pkin(2) * t10 + t29 * t39 + t36;
t31 = -t23 * t28 + t25 * t44;
t30 = -t28 * qJ(3) + t37;
t14 = qJ(6) + t21;
t13 = cos(t21);
t8 = cos(t14);
t7 = sin(t14);
t1 = pkin(5) * t13 + t11;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t26, 0, 0; t26, t29, 0, 0; 0, 0, 1, t22; 0, 0, 0, 1; t10, -t41, t26, t36; t9, -t43, -t29, t33; t25, t28, 0, t22; 0, 0, 0, 1; t26, -t10, t41, t32; -t29, -t9, t43, t33 + t45; 0, -t25, -t28, t30; 0, 0, 0, 1; t24 * t41 + t42, -t26 * t24 + t25 * t40, t10, t26 * pkin(3) + pkin(8) * t10 + t32; t24 * t43 - t40, t29 * t24 + t25 * t42, t9, pkin(8) * t9 + (-pkin(3) - pkin(7)) * t29 + t35; -t28 * t24, -t28 * t27, t25, t25 * pkin(8) + t30; 0, 0, 0, 1; t12 * t41 + t26 * t13, -t26 * t12 + t13 * t41, t10, t26 * t11 + t31 * t29 + t32; t12 * t43 - t29 * t13, t29 * t12 + t13 * t43, t9 (-pkin(7) - t11) * t29 + t31 * t26 + t35; -t28 * t12, -t28 * t13, t25, -t25 * t23 + (-qJ(3) - t44) * t28 + t37; 0, 0, 0, 1; t26 * t8 + t7 * t41, -t26 * t7 + t8 * t41, t10, t26 * t1 + t34 * t29 + t32; -t29 * t8 + t7 * t43, t29 * t7 + t8 * t43, t9 (-pkin(7) - t1) * t29 + t34 * t26 + t35; -t28 * t7, -t28 * t8, t25, -t25 * t20 + (-qJ(3) - t2) * t28 + t37; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
