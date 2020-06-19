% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% 
% Output:
% Tc_mdh [4x4x(2+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   3:  mdh base (link 0) -> mdh frame (3-1), link (3-1)
%   ...
%   2+1:  mdh base (link 0) -> mdh frame (2)
% T_c_stack [(2+1)*3 x 4]
%   stacked matrices from Tc_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S2RR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:17
% EndTime: 2020-06-19 09:14:17
% DurationCPUTime: 0.08s
% Computational Cost: add. (11->7), mult. (2->2), div. (0->0), fcn. (10->4), ass. (0->7)
t6 = pkin(2) + 0;
t5 = cos(qJ(1));
t4 = sin(qJ(1));
t3 = qJ(1) + qJ(2);
t2 = cos(t3);
t1 = sin(t3);
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t5, -t4, 0, 0; t4, t5, 0, 0; 0, 0, 1, t6; t2, -t1, 0, t5 * pkin(1) + 0; t1, t2, 0, t4 * pkin(1) + 0; 0, 0, 1, pkin(3) + t6;];
Tc_stack = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,2+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,2+1]); end % symbolisch
for i = 1:2+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
