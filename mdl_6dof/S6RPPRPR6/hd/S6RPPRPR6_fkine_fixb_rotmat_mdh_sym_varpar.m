% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% Datum: 2018-11-23 15:42
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPPRPR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:42:15
% EndTime: 2018-11-23 15:42:15
% DurationCPUTime: 0.10s
% Computational Cost: add. (81->46), mult. (85->35), div. (0->0), fcn. (126->6), ass. (0->28)
t15 = sin(qJ(4));
t16 = sin(qJ(1));
t3 = t16 * t15;
t18 = cos(qJ(4));
t33 = t16 * t18;
t19 = cos(qJ(1));
t4 = t19 * t15;
t32 = t19 * t18;
t31 = qJ(5) * t18;
t13 = pkin(6) + 0;
t30 = t16 * pkin(1) + 0;
t29 = pkin(2) + t13;
t28 = t19 * pkin(1) + t16 * qJ(2) + 0;
t27 = pkin(3) + t29;
t26 = t19 * qJ(3) + t28;
t25 = pkin(8) * t15 - t31;
t24 = -t19 * qJ(2) + t30;
t23 = t18 * pkin(4) + t15 * qJ(5) + t27;
t6 = t16 * qJ(3);
t22 = t24 + t6;
t21 = -t16 * pkin(7) + t26;
t11 = t19 * pkin(7);
t20 = t11 + t22;
t17 = cos(qJ(6));
t14 = sin(qJ(6));
t2 = pkin(4) * t4;
t1 = pkin(4) * t3;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t19, -t16, 0, 0; t16, t19, 0, 0; 0, 0, 1, t13; 0, 0, 0, 1; 0, -t19, t16, t28; 0, -t16, -t19, t24; 1, 0, 0, t13; 0, 0, 0, 1; 0, t16, t19, t26; 0, -t19, t16, t22; 1, 0, 0, t29; 0, 0, 0, 1; t4, t32, -t16, t21; t3, t33, t19, t20; t18, -t15, 0, t27; 0, 0, 0, 1; -t16, -t4, -t32, -t19 * t31 + t2 + t21; t19, -t3, -t33, -t16 * t31 + t1 + t20; 0, -t18, t15, t23; 0, 0, 0, 1; -t14 * t32 - t16 * t17, t16 * t14 - t17 * t32, t4, t2 + t25 * t19 + (-pkin(5) - pkin(7)) * t16 + t26; -t14 * t33 + t19 * t17, -t19 * t14 - t17 * t33, t3, t1 + t11 + t6 + (pkin(5) - qJ(2)) * t19 + t25 * t16 + t30; t15 * t14, t15 * t17, t18, t18 * pkin(8) + t23; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
