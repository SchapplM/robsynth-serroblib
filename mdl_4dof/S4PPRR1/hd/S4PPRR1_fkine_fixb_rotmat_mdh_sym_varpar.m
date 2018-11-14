% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% T_c_mdh [4x4x(4+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   4+1:  mdh base (link 0) -> mdh frame (4)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:40
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S4PPRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:06
% EndTime: 2018-11-14 13:40:06
% DurationCPUTime: 0.07s
% Computational Cost: add. (45->24), mult. (38->20), div. (0->0), fcn. (62->6), ass. (0->19)
t13 = sin(pkin(6));
t15 = sin(qJ(3));
t21 = t13 * t15;
t20 = t13 * pkin(1) + 0;
t11 = qJ(1) + 0;
t14 = cos(pkin(6));
t19 = t14 * pkin(1) + t13 * qJ(2) + 0;
t18 = -pkin(4) + t11;
t17 = -t14 * qJ(2) + t20;
t16 = cos(qJ(3));
t12 = qJ(3) + qJ(4);
t8 = cos(t12);
t7 = sin(t12);
t5 = t16 * pkin(3) + pkin(2);
t4 = t13 * t16 - t14 * t15;
t3 = -t14 * t16 - t21;
t2 = t13 * t8 - t14 * t7;
t1 = -t13 * t7 - t14 * t8;
t6 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t14, -t13, 0, 0; t13, t14, 0, 0; 0, 0, 1, t11; 0, 0, 0, 1; t14, 0, t13, t19; t13, 0, -t14, t17; 0, 1, 0, t11; 0, 0, 0, 1; -t3, t4, 0, t14 * pkin(2) + t19; t4, t3, 0, t13 * pkin(2) + t17; 0, 0, -1, t18; 0, 0, 0, 1; -t1, t2, 0, pkin(3) * t21 + t14 * t5 + t19; t2, t1, 0, t13 * t5 + (-pkin(3) * t15 - qJ(2)) * t14 + t20; 0, 0, -1, -pkin(5) + t18; 0, 0, 0, 1;];
T_ges = t6;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
