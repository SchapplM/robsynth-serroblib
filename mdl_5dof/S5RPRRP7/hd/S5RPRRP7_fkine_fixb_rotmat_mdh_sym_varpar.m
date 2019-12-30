% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-29 17:21
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRRP7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 17:20:27
% EndTime: 2019-12-29 17:20:27
% DurationCPUTime: 0.19s
% Computational Cost: add. (125->35), mult. (96->35), div. (0->0), fcn. (142->8), ass. (0->31)
t22 = qJ(1) + pkin(8);
t16 = sin(t22);
t24 = sin(qJ(3));
t9 = t16 * t24;
t27 = cos(qJ(3));
t41 = t16 * t27;
t17 = cos(t22);
t11 = t17 * t24;
t40 = t17 * t27;
t23 = sin(qJ(4));
t39 = t23 * t27;
t38 = t24 * t23;
t26 = cos(qJ(4));
t37 = t26 * t27;
t36 = pkin(5) + 0;
t25 = sin(qJ(1));
t35 = t25 * pkin(1) + 0;
t28 = cos(qJ(1));
t34 = t28 * pkin(1) + 0;
t18 = qJ(2) + t36;
t33 = t17 * pkin(2) + t16 * pkin(6) + t34;
t32 = t16 * pkin(2) - t17 * pkin(6) + t35;
t31 = pkin(3) * t40 + pkin(7) * t11 + t33;
t30 = t24 * pkin(3) - t27 * pkin(7) + t18;
t29 = pkin(3) * t41 + pkin(7) * t9 + t32;
t15 = t24 * t26;
t4 = t16 * t23 + t17 * t37;
t3 = -t16 * t26 + t17 * t39;
t2 = t16 * t37 - t17 * t23;
t1 = t16 * t39 + t17 * t26;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t25, 0, 0; t25, t28, 0, 0; 0, 0, 1, t36; 0, 0, 0, 1; t17, -t16, 0, t34; t16, t17, 0, t35; 0, 0, 1, t18; 0, 0, 0, 1; t40, -t11, t16, t33; t41, -t9, -t17, t32; t24, t27, 0, t18; 0, 0, 0, 1; t4, -t3, t11, t31; t2, -t1, t9, t29; t15, -t38, -t27, t30; 0, 0, 0, 1; t4, t11, t3, t4 * pkin(4) + t3 * qJ(5) + t31; t2, t9, t1, t2 * pkin(4) + t1 * qJ(5) + t29; t15, -t27, t38, (pkin(4) * t26 + qJ(5) * t23) * t24 + t30; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
