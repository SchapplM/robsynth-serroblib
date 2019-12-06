% Calculate homogenous joint transformation matrices for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S5PRRRR10_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:05
% EndTime: 2019-12-05 17:23:05
% DurationCPUTime: 0.04s
% Computational Cost: add. (11->11), mult. (12->12), div. (0->0), fcn. (40->14), ass. (0->15)
t90 = cos(qJ(2));
t89 = cos(qJ(3));
t88 = cos(qJ(4));
t87 = cos(qJ(5));
t86 = sin(qJ(2));
t85 = sin(qJ(3));
t84 = sin(qJ(4));
t83 = sin(qJ(5));
t82 = cos(pkin(5));
t81 = cos(pkin(6));
t80 = cos(pkin(11));
t79 = sin(pkin(5));
t78 = sin(pkin(6));
t77 = sin(pkin(11));
t1 = [t80, -t77, 0, 0; t77, t80, 0, 0; 0, 0, 1, qJ(1); 0, 0, 0, 1; t90, -t86, 0, pkin(1); t82 * t86, t82 * t90, -t79, -t79 * pkin(7); t79 * t86, t79 * t90, t82, t82 * pkin(7); 0, 0, 0, 1; t89, -t85, 0, pkin(2); t81 * t85, t81 * t89, -t78, -t78 * pkin(8); t78 * t85, t78 * t89, t81, t81 * pkin(8); 0, 0, 0, 1; t88, -t84, 0, pkin(3); 0, 0, -1, -pkin(9); t84, t88, 0, 0; 0, 0, 0, 1; t87, -t83, 0, pkin(4); 0, 0, -1, -pkin(10); t83, t87, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
