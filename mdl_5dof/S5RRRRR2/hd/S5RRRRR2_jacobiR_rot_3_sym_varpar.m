% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S5RRRRR2
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
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRRR2_jacobiR_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_jacobiR_rot_3_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_jacobiR_rot_3_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:26:51
% EndTime: 2019-03-29 15:26:51
% DurationCPUTime: 0.02s
% Computational Cost: add. (25->10), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
t21 = qJ(1) + qJ(2);
t19 = sin(t21);
t23 = cos(qJ(3));
t25 = t19 * t23;
t20 = cos(t21);
t22 = sin(qJ(3));
t24 = t20 * t22;
t18 = t20 * t23;
t17 = t19 * t22;
t1 = [-t25, -t25, -t24, 0, 0; t18, t18, -t17, 0, 0; 0, 0, t23, 0, 0; t17, t17, -t18, 0, 0; -t24, -t24, -t25, 0, 0; 0, 0, -t22, 0, 0; t20, t20, 0, 0, 0; t19, t19, 0, 0, 0; 0, 0, 0, 0, 0;];
JR_rot  = t1;
