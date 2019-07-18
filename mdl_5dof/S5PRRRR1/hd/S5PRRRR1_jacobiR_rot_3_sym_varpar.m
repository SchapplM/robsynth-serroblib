% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRRRR1_jacobiR_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobiR_rot_3_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobiR_rot_3_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:29:14
% EndTime: 2019-07-18 13:29:14
% DurationCPUTime: 0.02s
% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
t20 = sin(qJ(3));
t21 = sin(qJ(2));
t27 = t21 * t20;
t22 = cos(qJ(3));
t26 = t21 * t22;
t23 = cos(qJ(2));
t25 = t23 * t20;
t24 = t23 * t22;
t1 = [0, -t26, -t25, 0, 0; 0, 0, -t22, 0, 0; 0, t24, -t27, 0, 0; 0, t27, -t24, 0, 0; 0, 0, t20, 0, 0; 0, -t25, -t26, 0, 0; 0, t23, 0, 0, 0; 0, 0, 0, 0, 0; 0, t21, 0, 0, 0;];
JR_rot  = t1;
