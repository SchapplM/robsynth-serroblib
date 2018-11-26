% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2018-11-23 17:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPRRR11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:28:24
% EndTime: 2018-11-23 17:28:24
% DurationCPUTime: 0.16s
% Computational Cost: add. (151->70), mult. (162->70), div. (0->0), fcn. (221->10), ass. (0->37)
t25 = sin(qJ(1));
t24 = sin(qJ(2));
t39 = qJ(3) * t24;
t27 = cos(qJ(2));
t7 = t25 * t27;
t45 = pkin(2) * t7 + t25 * t39;
t29 = -pkin(9) - pkin(8);
t23 = sin(qJ(4));
t44 = t23 * pkin(4);
t26 = cos(qJ(4));
t11 = t26 * pkin(4) + pkin(3);
t43 = t25 * t24;
t42 = t25 * t26;
t28 = cos(qJ(1));
t41 = t28 * t24;
t40 = t28 * t26;
t8 = t28 * t27;
t22 = qJ(4) + qJ(5);
t20 = pkin(6) + 0;
t38 = t25 * pkin(1) + 0;
t37 = t24 * pkin(2) + t20;
t36 = t28 * pkin(1) + t25 * pkin(7) + 0;
t35 = t38 + t45;
t12 = sin(t22);
t2 = pkin(5) * t12 + t44;
t21 = -pkin(10) + t29;
t34 = t2 * t24 - t21 * t27;
t33 = -t28 * pkin(7) + t38;
t32 = pkin(2) * t8 + t28 * t39 + t36;
t31 = t24 * t44 - t27 * t29;
t30 = -t27 * qJ(3) + t37;
t14 = qJ(6) + t22;
t13 = cos(t22);
t10 = cos(t14);
t9 = sin(t14);
t1 = pkin(5) * t13 + t11;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t25, 0, 0; t25, t28, 0, 0; 0, 0, 1, t20; 0, 0, 0, 1; t8, -t41, t25, t36; t7, -t43, -t28, t33; t24, t27, 0, t20; 0, 0, 0, 1; t25, -t8, t41, t32; -t28, -t7, t43, t33 + t45; 0, -t24, -t27, t30; 0, 0, 0, 1; t23 * t41 + t42, -t25 * t23 + t24 * t40, t8, t25 * pkin(3) + pkin(8) * t8 + t32; t23 * t43 - t40, t28 * t23 + t24 * t42, t7, pkin(8) * t7 + (-pkin(3) - pkin(7)) * t28 + t35; -t27 * t23, -t27 * t26, t24, t24 * pkin(8) + t30; 0, 0, 0, 1; t12 * t41 + t25 * t13, -t25 * t12 + t13 * t41, t8, t25 * t11 + t31 * t28 + t32; t12 * t43 - t28 * t13, t28 * t12 + t13 * t43, t7 (-pkin(7) - t11) * t28 + t31 * t25 + t35; -t27 * t12, -t27 * t13, t24, -t24 * t29 + (-qJ(3) - t44) * t27 + t37; 0, 0, 0, 1; t25 * t10 + t9 * t41, t10 * t41 - t25 * t9, t8, t25 * t1 + t34 * t28 + t32; -t28 * t10 + t9 * t43, t10 * t43 + t28 * t9, t7 (-pkin(7) - t1) * t28 + t34 * t25 + t35; -t27 * t9, -t27 * t10, t24, -t24 * t21 + (-qJ(3) - t2) * t27 + t37; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
