% Calculate homogenous joint transformation matrices for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2022-02-03 23:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = S5RRRPR10_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-02-03 23:44:03
% EndTime: 2022-02-03 23:44:03
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (6->6), div. (0->0), fcn. (30->12), ass. (0->13)
t78 = cos(qJ(1));
t77 = cos(qJ(2));
t76 = cos(qJ(3));
t75 = cos(qJ(5));
t74 = sin(qJ(1));
t73 = sin(qJ(2));
t72 = sin(qJ(3));
t71 = sin(qJ(5));
t70 = cos(pkin(5));
t69 = cos(pkin(10));
t68 = sin(pkin(5));
t67 = sin(pkin(10));
t1 = [t78, -t74, 0, 0; t74, t78, 0, 0; 0, 0, 1, pkin(6); t77, -t73, 0, pkin(1); t70 * t73, t70 * t77, -t68, -t68 * pkin(7); t68 * t73, t68 * t77, t70, t70 * pkin(7); t76, -t72, 0, pkin(2); 0, 0, -1, -pkin(8); t72, t76, 0, 0; t69, -t67, 0, pkin(3); t67, t69, 0, 0; 0, 0, 1, qJ(4); t75, -t71, 0, pkin(4); 0, 0, -1, -pkin(9); t71, t75, 0, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
