% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2018-11-23 16:19
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRPR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:19:37
% EndTime: 2018-11-23 16:19:37
% DurationCPUTime: 0.11s
% Computational Cost: add. (142->50), mult. (87->46), div. (0->0), fcn. (132->10), ass. (0->35)
t22 = -pkin(8) - pkin(7);
t17 = sin(qJ(3));
t40 = t17 * pkin(3);
t18 = sin(qJ(1));
t15 = qJ(3) + qJ(4);
t6 = pkin(10) + t15;
t4 = cos(t6);
t39 = t18 * t4;
t21 = cos(qJ(1));
t38 = t21 * t4;
t16 = sin(qJ(6));
t37 = t18 * t16;
t36 = t18 * t17;
t19 = cos(qJ(6));
t35 = t18 * t19;
t34 = t21 * t16;
t33 = t21 * t19;
t7 = sin(t15);
t2 = pkin(4) * t7 + t40;
t32 = -qJ(2) - t2;
t14 = pkin(6) + 0;
t31 = t18 * pkin(1) + 0;
t30 = pkin(2) + t14;
t29 = t21 * pkin(1) + t18 * qJ(2) + 0;
t20 = cos(qJ(3));
t28 = t20 * pkin(3) + t30;
t3 = sin(t6);
t27 = pkin(5) * t3 - pkin(9) * t4;
t8 = cos(t15);
t26 = pkin(4) * t8 + t28;
t13 = -qJ(5) + t22;
t25 = -t18 * t13 + t31;
t24 = -t21 * qJ(2) + t31;
t23 = -t21 * t13 + t18 * t2 + t29;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t21, -t18, 0, 0; t18, t21, 0, 0; 0, 0, 1, t14; 0, 0, 0, 1; 0, -t21, t18, t29; 0, -t18, -t21, t24; 1, 0, 0, t14; 0, 0, 0, 1; t36, t18 * t20, t21, t21 * pkin(7) + t29; -t21 * t17, -t21 * t20, t18, t18 * pkin(7) + t24; t20, -t17, 0, t30; 0, 0, 0, 1; t18 * t7, t18 * t8, t21, pkin(3) * t36 - t21 * t22 + t29; -t21 * t7, -t21 * t8, t18, -t18 * t22 + (-qJ(2) - t40) * t21 + t31; t8, -t7, 0, t28; 0, 0, 0, 1; t18 * t3, t39, t21, t23; -t21 * t3, -t38, t18, t32 * t21 + t25; t4, -t3, 0, t26; 0, 0, 0, 1; t3 * t35 + t34, -t3 * t37 + t33, -t39, t27 * t18 + t23; -t3 * t33 + t37, t3 * t34 + t35, t38 (-t27 + t32) * t21 + t25; t4 * t19, -t4 * t16, t3, t4 * pkin(5) + t3 * pkin(9) + t26; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
