% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2018-11-23 16:46
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPPRP3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:45:57
% EndTime: 2018-11-23 16:45:57
% DurationCPUTime: 0.12s
% Computational Cost: add. (106->54), mult. (152->46), div. (0->0), fcn. (206->6), ass. (0->40)
t27 = sin(qJ(1));
t29 = cos(qJ(2));
t12 = t27 * t29;
t26 = sin(qJ(2));
t44 = qJ(3) * t26;
t50 = pkin(2) * t12 + t27 * t44;
t30 = cos(qJ(1));
t49 = pkin(3) * t12 + t30 * qJ(4);
t25 = sin(qJ(5));
t48 = pkin(5) * t25;
t11 = t27 * t26;
t28 = cos(qJ(5));
t47 = t27 * t28;
t46 = t29 * t28;
t14 = t30 * t26;
t45 = t30 * t28;
t15 = t30 * t29;
t23 = pkin(6) + 0;
t43 = t27 * pkin(1) + 0;
t42 = t26 * pkin(2) + t23;
t41 = t30 * pkin(1) + t27 * pkin(7) + 0;
t18 = t26 * pkin(3);
t40 = t18 + t42;
t39 = pkin(4) * t26 + pkin(8) * t29;
t16 = t28 * pkin(5) + pkin(4);
t24 = -qJ(6) - pkin(8);
t38 = t16 * t26 - t24 * t29;
t37 = -t30 * pkin(7) + t43;
t36 = pkin(2) * t15 + t30 * t44 + t41;
t35 = pkin(3) * t15 + t36;
t34 = -t29 * qJ(3) + t42;
t33 = t37 + t50;
t32 = t33 + t49;
t31 = -t27 * qJ(4) + t35;
t13 = t29 * t25;
t4 = -t27 * t25 + t26 * t45;
t3 = -t25 * t14 - t47;
t2 = t30 * t25 + t26 * t47;
t1 = -t25 * t11 + t45;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t30, -t27, 0, 0; t27, t30, 0, 0; 0, 0, 1, t23; 0, 0, 0, 1; t15, -t14, t27, t41; t12, -t11, -t30, t37; t26, t29, 0, t23; 0, 0, 0, 1; t15, t27, t14, t36; t12, -t30, t11, t33; t26, 0, -t29, t34; 0, 0, 0, 1; t14, -t15, -t27, t31; t11, -t12, t30, t32; -t29, -t26, 0, t18 + t34; 0, 0, 0, 1; t4, t3, t15, t39 * t30 + t31; t2, t1, t12, t39 * t27 + t32; -t46, t13, t26, t26 * pkin(8) + (-pkin(4) - qJ(3)) * t29 + t40; 0, 0, 0, 1; t4, t3, t15, t38 * t30 + (-qJ(4) - t48) * t27 + t35; t2, t1, t12 (-pkin(7) + t48) * t30 + t38 * t27 + t43 + t49 + t50; -t46, t13, t26, -t26 * t24 + (-qJ(3) - t16) * t29 + t40; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
