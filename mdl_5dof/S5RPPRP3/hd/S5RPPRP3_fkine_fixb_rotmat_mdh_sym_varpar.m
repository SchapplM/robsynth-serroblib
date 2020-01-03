% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPPRP3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:38
% EndTime: 2019-12-31 17:50:39
% DurationCPUTime: 0.08s
% Computational Cost: add. (89->30), mult. (37->18), div. (0->0), fcn. (65->6), ass. (0->21)
t13 = sin(qJ(4));
t11 = qJ(1) + pkin(7);
t6 = sin(t11);
t1 = t6 * t13;
t7 = cos(t11);
t25 = t7 * t13;
t15 = cos(qJ(4));
t24 = t7 * t15;
t23 = pkin(5) + 0;
t14 = sin(qJ(1));
t22 = t14 * pkin(1) + 0;
t16 = cos(qJ(1));
t21 = t16 * pkin(1) + 0;
t20 = t6 * pkin(2) + t22;
t8 = qJ(2) + t23;
t19 = t7 * pkin(2) + t6 * qJ(3) + t21;
t18 = pkin(3) + t8;
t17 = -t7 * qJ(3) + t20;
t12 = -qJ(5) - pkin(6);
t2 = t6 * t15;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t16, -t14, 0, 0; t14, t16, 0, 0; 0, 0, 1, t23; 0, 0, 0, 1; t7, -t6, 0, t21; t6, t7, 0, t22; 0, 0, 1, t8; 0, 0, 0, 1; 0, -t7, t6, t19; 0, -t6, -t7, t17; 1, 0, 0, t8; 0, 0, 0, 1; t1, t2, t7, t7 * pkin(6) + t19; -t25, -t24, t6, t6 * pkin(6) + t17; t15, -t13, 0, t18; 0, 0, 0, 1; t1, t2, t7, pkin(4) * t1 - t7 * t12 + t19; -t25, -t24, t6, -t6 * t12 + (-pkin(4) * t13 - qJ(3)) * t7 + t20; t15, -t13, 0, t15 * pkin(4) + t18; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
