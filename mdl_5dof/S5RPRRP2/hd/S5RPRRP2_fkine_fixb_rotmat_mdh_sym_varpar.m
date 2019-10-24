% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP2
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
% Datum: 2019-10-24 10:43
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRRP2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:43:39
% EndTime: 2019-10-24 10:43:39
% DurationCPUTime: 0.08s
% Computational Cost: add. (110->33), mult. (33->18), div. (0->0), fcn. (61->8), ass. (0->24)
t16 = cos(qJ(4));
t12 = qJ(1) + pkin(8);
t10 = qJ(3) + t12;
t4 = sin(t10);
t25 = t4 * t16;
t14 = sin(qJ(4));
t5 = cos(t10);
t24 = t5 * t14;
t23 = pkin(5) + 0;
t17 = cos(qJ(1));
t22 = t17 * pkin(1) + 0;
t9 = cos(t12);
t21 = pkin(2) * t9 + t22;
t20 = qJ(2) + t23;
t15 = sin(qJ(1));
t19 = -t15 * pkin(1) + 0;
t7 = pkin(6) + t20;
t8 = sin(t12);
t18 = -pkin(2) * t8 + t19;
t13 = -qJ(5) - pkin(7);
t6 = t16 * pkin(4) + pkin(3);
t2 = t5 * t16;
t1 = t4 * t14;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t23; -t15, -t17, 0, 0; t17, -t15, 0, 0; 0, 0, 0, 1; 0, 0, 1, t20; -t8, -t9, 0, t19; t9, -t8, 0, t22; 0, 0, 0, 1; 0, 0, 1, t7; -t4, -t5, 0, t18; t5, -t4, 0, t21; 0, 0, 0, 1; t14, t16, 0, t7; -t25, t1, t5, -t4 * pkin(3) + t5 * pkin(7) + t18; t2, -t24, t4, t5 * pkin(3) + t4 * pkin(7) + t21; 0, 0, 0, 1; t14, t16, 0, t14 * pkin(4) + t7; -t25, t1, t5, -t5 * t13 - t4 * t6 + t18; t2, -t24, t4, -t4 * t13 + t5 * t6 + t21; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
