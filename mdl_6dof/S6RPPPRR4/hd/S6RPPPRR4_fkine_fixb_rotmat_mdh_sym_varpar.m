% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% Datum: 2018-11-23 15:38
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPPPRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:38:29
% EndTime: 2018-11-23 15:38:29
% DurationCPUTime: 0.10s
% Computational Cost: add. (117->43), mult. (162->38), div. (0->0), fcn. (243->8), ass. (0->29)
t29 = sin(pkin(9));
t30 = cos(pkin(9));
t33 = sin(qJ(1));
t34 = cos(qJ(1));
t3 = -t33 * t29 - t34 * t30;
t38 = t3 * pkin(7);
t4 = t34 * t29 - t33 * t30;
t37 = t4 * pkin(7);
t18 = cos(qJ(5));
t36 = t3 * t18;
t35 = t4 * t18;
t15 = sin(qJ(6));
t16 = sin(qJ(5));
t32 = t15 * t16;
t17 = cos(qJ(6));
t31 = t16 * t17;
t14 = pkin(6) + 0;
t28 = t34 * pkin(1) + t33 * qJ(2) + 0;
t8 = -qJ(3) + t14;
t27 = -pkin(4) + t8;
t26 = t34 * pkin(2) + t28;
t25 = -t3 * pkin(3) + t26;
t24 = pkin(5) * t16 - pkin(8) * t18 + qJ(4);
t23 = t33 * pkin(1) - t34 * qJ(2) + 0;
t22 = t33 * pkin(2) + t23;
t21 = t4 * qJ(4) + t25;
t20 = -t4 * pkin(3) + t22;
t19 = -t3 * qJ(4) + t20;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t33, 0, 0; t33, t34, 0, 0; 0, 0, 1, t14; 0, 0, 0, 1; t34, 0, t33, t28; t33, 0, -t34, t23; 0, 1, 0, t14; 0, 0, 0, 1; -t3, -t4, 0, t26; -t4, t3, 0, t22; 0, 0, -1, t8; 0, 0, 0, 1; 0, t3, t4, t21; 0, t4, -t3, t19; -1, 0, 0, t8; 0, 0, 0, 1; t4 * t16, t35, -t3, t21 - t38; -t3 * t16, -t36, -t4, t19 - t37; -t18, t16, 0, t27; 0, 0, 0, 1; -t3 * t15 + t4 * t31, -t3 * t17 - t4 * t32, -t35, t24 * t4 + t25 - t38; -t4 * t15 - t3 * t31, -t4 * t17 + t3 * t32, t36, -t24 * t3 + t20 - t37; -t18 * t17, t18 * t15, -t16, -t18 * pkin(5) - t16 * pkin(8) + t27; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
