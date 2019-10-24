% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-10-24 10:32
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRRRP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:32:36
% EndTime: 2019-10-24 10:32:36
% DurationCPUTime: 0.08s
% Computational Cost: add. (106->29), mult. (33->18), div. (0->0), fcn. (61->8), ass. (0->24)
t18 = sin(qJ(4));
t14 = pkin(8) + qJ(2);
t11 = qJ(3) + t14;
t5 = sin(t11);
t27 = t5 * t18;
t6 = cos(t11);
t26 = t6 * t18;
t15 = sin(pkin(8));
t25 = t15 * pkin(1) + 0;
t16 = cos(pkin(8));
t24 = t16 * pkin(1) + 0;
t23 = qJ(1) + 0;
t9 = sin(t14);
t22 = pkin(2) * t9 + t25;
t10 = cos(t14);
t21 = pkin(2) * t10 + t24;
t20 = pkin(5) + t23;
t8 = pkin(6) + t20;
t19 = cos(qJ(4));
t17 = -qJ(5) - pkin(7);
t7 = t19 * pkin(4) + pkin(3);
t2 = t6 * t19;
t1 = t5 * t19;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t16, -t15, 0, 0; t15, t16, 0, 0; 0, 0, 1, t23; 0, 0, 0, 1; t10, -t9, 0, t24; t9, t10, 0, t25; 0, 0, 1, t20; 0, 0, 0, 1; t6, -t5, 0, t21; t5, t6, 0, t22; 0, 0, 1, t8; 0, 0, 0, 1; t2, -t26, t5, t6 * pkin(3) + t5 * pkin(7) + t21; t1, -t27, -t6, t5 * pkin(3) - t6 * pkin(7) + t22; t18, t19, 0, t8; 0, 0, 0, 1; t2, -t26, t5, -t5 * t17 + t6 * t7 + t21; t1, -t27, -t6, t6 * t17 + t5 * t7 + t22; t18, t19, 0, t18 * pkin(4) + t8; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
