% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPRR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:23
% EndTime: 2019-12-31 20:15:23
% DurationCPUTime: 0.08s
% Computational Cost: add. (95->31), mult. (37->22), div. (0->0), fcn. (65->8), ass. (0->21)
t13 = sin(qJ(4));
t12 = qJ(1) + qJ(2);
t5 = sin(t12);
t25 = t5 * t13;
t24 = pkin(5) + 0;
t14 = sin(qJ(1));
t23 = t14 * pkin(1) + 0;
t16 = cos(qJ(1));
t22 = t16 * pkin(1) + 0;
t8 = pkin(6) + t24;
t21 = t5 * pkin(2) + t23;
t20 = pkin(3) + t8;
t7 = cos(t12);
t19 = t7 * pkin(2) + t5 * qJ(3) + t22;
t18 = -t7 * qJ(3) + t21;
t17 = -pkin(8) - pkin(7);
t15 = cos(qJ(4));
t11 = qJ(4) + qJ(5);
t6 = cos(t11);
t4 = sin(t11);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t16, -t14, 0, 0; t14, t16, 0, 0; 0, 0, 1, t24; 0, 0, 0, 1; t7, -t5, 0, t22; t5, t7, 0, t23; 0, 0, 1, t8; 0, 0, 0, 1; 0, -t7, t5, t19; 0, -t5, -t7, t18; 1, 0, 0, t8; 0, 0, 0, 1; t25, t5 * t15, t7, t7 * pkin(7) + t19; -t7 * t13, -t7 * t15, t5, t5 * pkin(7) + t18; t15, -t13, 0, t20; 0, 0, 0, 1; t5 * t4, t5 * t6, t7, pkin(4) * t25 - t7 * t17 + t19; -t7 * t4, -t7 * t6, t5, -t5 * t17 + (-pkin(4) * t13 - qJ(3)) * t7 + t21; t6, -t4, 0, t15 * pkin(4) + t20; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
