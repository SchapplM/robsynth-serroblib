% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRPP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:02
% EndTime: 2019-12-31 18:16:02
% DurationCPUTime: 0.08s
% Computational Cost: add. (60->37), mult. (61->25), div. (0->0), fcn. (93->4), ass. (0->20)
t13 = sin(qJ(3));
t14 = sin(qJ(1));
t3 = t14 * t13;
t15 = cos(qJ(3));
t27 = t14 * t15;
t16 = cos(qJ(1));
t26 = t16 * t13;
t25 = qJ(4) * t15;
t12 = pkin(5) + 0;
t24 = t14 * pkin(1) + 0;
t23 = pkin(2) + t12;
t22 = t16 * pkin(1) + t14 * qJ(2) + 0;
t7 = t14 * pkin(6);
t21 = t16 * t25 + t24 + t7;
t20 = t16 * pkin(6) + t22;
t19 = t15 * pkin(3) + t13 * qJ(4) + t23;
t18 = pkin(3) * t3 + t20;
t17 = -t16 * qJ(2) + t24;
t4 = t16 * t15;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t16, -t14, 0, 0; t14, t16, 0, 0; 0, 0, 1, t12; 0, 0, 0, 1; 0, -t16, t14, t22; 0, -t14, -t16, t17; 1, 0, 0, t12; 0, 0, 0, 1; t3, t27, t16, t20; -t26, -t4, t14, t17 + t7; t15, -t13, 0, t23; 0, 0, 0, 1; t3, t16, -t27, -t14 * t25 + t18; -t26, t14, t4, (-pkin(3) * t13 - qJ(2)) * t16 + t21; t15, 0, t13, t19; 0, 0, 0, 1; t3, -t27, -t16, -t16 * qJ(5) + (pkin(4) * t13 - t25) * t14 + t18; -t26, t4, -t14, -t14 * qJ(5) + (-qJ(2) + (-pkin(3) - pkin(4)) * t13) * t16 + t21; t15, t13, 0, t15 * pkin(4) + t19; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
