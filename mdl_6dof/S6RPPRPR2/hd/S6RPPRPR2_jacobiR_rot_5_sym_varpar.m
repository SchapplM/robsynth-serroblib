% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:09
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRPR2_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:09:18
% EndTime: 2019-02-22 10:09:18
% DurationCPUTime: 0.02s
% Computational Cost: add. (23->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->11)
t52 = pkin(10) + qJ(4);
t48 = sin(t52);
t53 = qJ(1) + pkin(9);
t49 = sin(t53);
t55 = t49 * t48;
t50 = cos(t52);
t51 = cos(t53);
t54 = t51 * t50;
t47 = t51 * t48;
t46 = t49 * t50;
t1 = [t51, 0, 0, 0, 0, 0; t49, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t46, 0, 0, t47, 0, 0; -t54, 0, 0, t55, 0, 0; 0, 0, 0, -t50, 0, 0; -t55, 0, 0, t54, 0, 0; t47, 0, 0, t46, 0, 0; 0, 0, 0, t48, 0, 0;];
JR_rot  = t1;
