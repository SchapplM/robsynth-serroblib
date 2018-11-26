% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2018-11-23 16:45
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPPRP2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:45:35
% EndTime: 2018-11-23 16:45:35
% DurationCPUTime: 0.12s
% Computational Cost: add. (163->57), mult. (130->50), div. (0->0), fcn. (184->8), ass. (0->41)
t20 = qJ(2) + pkin(9);
t18 = cos(t20);
t29 = cos(qJ(1));
t12 = t29 * t18;
t17 = sin(t20);
t39 = qJ(4) * t17;
t49 = pkin(3) * t12 + t29 * t39;
t24 = sin(qJ(5));
t48 = pkin(5) * t24;
t47 = t18 * t24;
t27 = cos(qJ(5));
t46 = t18 * t27;
t26 = sin(qJ(1));
t45 = t26 * t17;
t11 = t26 * t18;
t44 = t26 * t24;
t43 = t26 * t27;
t42 = t29 * t17;
t41 = t29 * t24;
t40 = t29 * t27;
t21 = pkin(6) + 0;
t28 = cos(qJ(2));
t16 = t28 * pkin(2) + pkin(1);
t38 = t29 * t16 + 0;
t25 = sin(qJ(2));
t37 = t25 * pkin(2) + t21;
t23 = -qJ(3) - pkin(7);
t36 = t26 * t16 + t29 * t23 + 0;
t35 = t38 + t49;
t34 = t17 * pkin(3) + t37;
t33 = pkin(3) * t11 + t26 * t39 + t36;
t32 = -t26 * t23 + t38;
t22 = -qJ(6) - pkin(8);
t31 = t17 * t48 - t18 * t22;
t30 = -t18 * qJ(4) + t34;
t15 = t27 * pkin(5) + pkin(4);
t4 = t17 * t44 - t40;
t3 = t17 * t43 + t41;
t2 = t17 * t41 + t43;
t1 = t17 * t40 - t44;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t26, 0, 0; t26, t29, 0, 0; 0, 0, 1, t21; 0, 0, 0, 1; t29 * t28, -t29 * t25, t26, t29 * pkin(1) + t26 * pkin(7) + 0; t26 * t28, -t26 * t25, -t29, t26 * pkin(1) - t29 * pkin(7) + 0; t25, t28, 0, t21; 0, 0, 0, 1; t12, -t42, t26, t32; t11, -t45, -t29, t36; t17, t18, 0, t37; 0, 0, 0, 1; t26, -t12, t42, t32 + t49; -t29, -t11, t45, t33; 0, -t17, -t18, t30; 0, 0, 0, 1; t2, t1, t12, pkin(8) * t12 + (pkin(4) - t23) * t26 + t35; t4, t3, t11, -t29 * pkin(4) + pkin(8) * t11 + t33; -t47, -t46, t17, t17 * pkin(8) + t30; 0, 0, 0, 1; t2, t1, t12, t31 * t29 + (t15 - t23) * t26 + t35; t4, t3, t11, -t29 * t15 + t31 * t26 + t33; -t47, -t46, t17, -t17 * t22 + (-qJ(4) - t48) * t18 + t34; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
