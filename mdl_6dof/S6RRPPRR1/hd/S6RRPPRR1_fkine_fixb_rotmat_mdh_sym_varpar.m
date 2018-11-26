% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2018-11-23 16:48
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPPRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:48:15
% EndTime: 2018-11-23 16:48:15
% DurationCPUTime: 0.15s
% Computational Cost: add. (194->54), mult. (174->56), div. (0->0), fcn. (244->10), ass. (0->36)
t27 = qJ(2) + pkin(10);
t24 = cos(t27);
t37 = cos(qJ(1));
t18 = t37 * t24;
t23 = sin(t27);
t47 = qJ(4) * t23;
t50 = pkin(3) * t18 + t37 * t47;
t33 = sin(qJ(1));
t49 = t33 * t23;
t17 = t33 * t24;
t48 = t37 * t23;
t28 = pkin(6) + 0;
t36 = cos(qJ(2));
t22 = t36 * pkin(2) + pkin(1);
t46 = t37 * t22 + 0;
t32 = sin(qJ(2));
t45 = t32 * pkin(2) + t28;
t29 = -qJ(3) - pkin(7);
t44 = t33 * t22 + t37 * t29 + 0;
t31 = sin(qJ(5));
t35 = cos(qJ(5));
t5 = t23 * t31 + t24 * t35;
t43 = -t33 * t29 + t46;
t42 = pkin(3) * t17 + t33 * t47 + t44;
t41 = t23 * pkin(3) - t24 * qJ(4) + t45;
t40 = pkin(4) * t17 + t37 * pkin(8) + t42;
t39 = t23 * pkin(4) + t41;
t38 = pkin(4) * t18 + (-pkin(8) - t29) * t33 + t46 + t50;
t34 = cos(qJ(6));
t30 = sin(qJ(6));
t6 = t23 * t35 - t24 * t31;
t4 = t5 * t37;
t3 = t18 * t31 - t35 * t48;
t2 = t5 * t33;
t1 = t17 * t31 - t35 * t49;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t37, -t33, 0, 0; t33, t37, 0, 0; 0, 0, 1, t28; 0, 0, 0, 1; t37 * t36, -t37 * t32, t33, t37 * pkin(1) + t33 * pkin(7) + 0; t33 * t36, -t33 * t32, -t37, t33 * pkin(1) - t37 * pkin(7) + 0; t32, t36, 0, t28; 0, 0, 0, 1; t18, -t48, t33, t43; t17, -t49, -t37, t44; t23, t24, 0, t45; 0, 0, 0, 1; t18, t33, t48, t43 + t50; t17, -t37, t49, t42; t23, 0, -t24, t41; 0, 0, 0, 1; t4, -t3, -t33, t38; t2, -t1, t37, t40; t6, -t5, 0, t39; 0, 0, 0, 1; -t33 * t30 + t4 * t34, -t4 * t30 - t33 * t34, t3, t4 * pkin(5) + t3 * pkin(9) + t38; t2 * t34 + t37 * t30, -t2 * t30 + t37 * t34, t1, t2 * pkin(5) + t1 * pkin(9) + t40; t6 * t34, -t6 * t30, t5, t6 * pkin(5) + t5 * pkin(9) + t39; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
