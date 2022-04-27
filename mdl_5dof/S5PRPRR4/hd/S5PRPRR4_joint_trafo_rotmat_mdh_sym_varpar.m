% Calculate homogenous joint transformation matrices for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2022-02-01 01:16
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = S5PRPRR4_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-02-01 01:15:56
% EndTime: 2022-02-01 01:15:56
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (6->6), div. (0->0), fcn. (30->12), ass. (0->13)
t87 = cos(qJ(2));
t86 = cos(qJ(4));
t85 = cos(qJ(5));
t84 = sin(qJ(2));
t83 = sin(qJ(4));
t82 = sin(qJ(5));
t81 = cos(pkin(5));
t80 = cos(pkin(9));
t79 = cos(pkin(10));
t78 = sin(pkin(5));
t77 = sin(pkin(9));
t76 = sin(pkin(10));
t1 = [t80, -t77, 0, 0; t77, t80, 0, 0; 0, 0, 1, qJ(1); t87, -t84, 0, pkin(1); t81 * t84, t81 * t87, -t78, -t78 * pkin(6); t78 * t84, t78 * t87, t81, t81 * pkin(6); t79, -t76, 0, pkin(2); t76, t79, 0, 0; 0, 0, 1, qJ(3); t86, -t83, 0, pkin(3); 0, 0, -1, -pkin(7); t83, t86, 0, 0; t85, -t82, 0, pkin(4); 0, 0, -1, -pkin(8); t82, t85, 0, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
