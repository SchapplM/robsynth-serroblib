% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2018-11-23 16:06
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRPRR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:06:29
% EndTime: 2018-11-23 16:06:29
% DurationCPUTime: 0.13s
% Computational Cost: add. (142->50), mult. (87->46), div. (0->0), fcn. (132->10), ass. (0->35)
t18 = sin(qJ(3));
t40 = t18 * pkin(3);
t19 = sin(qJ(1));
t14 = qJ(3) + pkin(10);
t8 = qJ(5) + t14;
t5 = cos(t8);
t39 = t19 * t5;
t22 = cos(qJ(1));
t38 = t22 * t5;
t17 = sin(qJ(6));
t37 = t19 * t17;
t36 = t19 * t18;
t20 = cos(qJ(6));
t35 = t19 * t20;
t34 = t22 * t17;
t33 = t22 * t20;
t6 = sin(t14);
t2 = pkin(4) * t6 + t40;
t32 = -qJ(2) - t2;
t16 = -qJ(4) - pkin(7);
t15 = pkin(6) + 0;
t31 = t19 * pkin(1) + 0;
t30 = pkin(2) + t15;
t29 = t22 * pkin(1) + t19 * qJ(2) + 0;
t21 = cos(qJ(3));
t28 = t21 * pkin(3) + t30;
t4 = sin(t8);
t27 = pkin(5) * t4 - pkin(9) * t5;
t7 = cos(t14);
t26 = pkin(4) * t7 + t28;
t13 = -pkin(8) + t16;
t25 = -t19 * t13 + t31;
t24 = -t22 * qJ(2) + t31;
t23 = -t22 * t13 + t19 * t2 + t29;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t22, -t19, 0, 0; t19, t22, 0, 0; 0, 0, 1, t15; 0, 0, 0, 1; 0, -t22, t19, t29; 0, -t19, -t22, t24; 1, 0, 0, t15; 0, 0, 0, 1; t36, t19 * t21, t22, t22 * pkin(7) + t29; -t22 * t18, -t22 * t21, t19, t19 * pkin(7) + t24; t21, -t18, 0, t30; 0, 0, 0, 1; t19 * t6, t19 * t7, t22, pkin(3) * t36 - t22 * t16 + t29; -t22 * t6, -t22 * t7, t19, -t19 * t16 + (-qJ(2) - t40) * t22 + t31; t7, -t6, 0, t28; 0, 0, 0, 1; t19 * t4, t39, t22, t23; -t22 * t4, -t38, t19, t32 * t22 + t25; t5, -t4, 0, t26; 0, 0, 0, 1; t4 * t35 + t34, -t4 * t37 + t33, -t39, t27 * t19 + t23; -t4 * t33 + t37, t4 * t34 + t35, t38 (-t27 + t32) * t22 + t25; t5 * t20, -t5 * t17, t4, t5 * pkin(5) + t4 * pkin(9) + t26; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
