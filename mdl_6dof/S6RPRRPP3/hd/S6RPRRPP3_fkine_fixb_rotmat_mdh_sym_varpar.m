% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2018-11-23 16:12
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRPP3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:12:20
% EndTime: 2018-11-23 16:12:20
% DurationCPUTime: 0.13s
% Computational Cost: add. (195->47), mult. (163->40), div. (0->0), fcn. (228->8), ass. (0->35)
t29 = sin(qJ(4));
t30 = sin(qJ(3));
t20 = t30 * t29;
t32 = cos(qJ(4));
t21 = t30 * t32;
t50 = pkin(4) * t21 + qJ(5) * t20;
t28 = qJ(1) + pkin(9);
t22 = sin(t28);
t11 = t22 * t30;
t33 = cos(qJ(3));
t49 = t22 * t33;
t23 = cos(t28);
t14 = t23 * t30;
t48 = t23 * t33;
t47 = t29 * t33;
t46 = t32 * t33;
t45 = pkin(6) + 0;
t31 = sin(qJ(1));
t44 = t31 * pkin(1) + 0;
t34 = cos(qJ(1));
t43 = t34 * pkin(1) + 0;
t24 = qJ(2) + t45;
t42 = t30 * pkin(3) + t24;
t41 = t23 * pkin(2) + t22 * pkin(7) + t43;
t40 = t22 * pkin(2) - t23 * pkin(7) + t44;
t39 = pkin(3) * t48 + pkin(8) * t14 + t41;
t38 = -t33 * pkin(8) + t42;
t37 = pkin(3) * t49 + pkin(8) * t11 + t40;
t5 = -t22 * t32 + t23 * t47;
t6 = t22 * t29 + t23 * t46;
t36 = t6 * pkin(4) + t5 * qJ(5) + t39;
t3 = t22 * t47 + t23 * t32;
t4 = t22 * t46 - t23 * t29;
t35 = t4 * pkin(4) + t3 * qJ(5) + t37;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t31, 0, 0; t31, t34, 0, 0; 0, 0, 1, t45; 0, 0, 0, 1; t23, -t22, 0, t43; t22, t23, 0, t44; 0, 0, 1, t24; 0, 0, 0, 1; t48, -t14, t22, t41; t49, -t11, -t23, t40; t30, t33, 0, t24; 0, 0, 0, 1; t6, -t5, t14, t39; t4, -t3, t11, t37; t21, -t20, -t33, t38; 0, 0, 0, 1; t14, -t6, t5, t36; t11, -t4, t3, t35; -t33, -t21, t20, t38 + t50; 0, 0, 0, 1; t14, t5, t6, pkin(5) * t14 + t6 * qJ(6) + t36; t11, t3, t4, pkin(5) * t11 + t4 * qJ(6) + t35; -t33, t20, t21, qJ(6) * t21 + (-pkin(5) - pkin(8)) * t33 + t42 + t50; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
