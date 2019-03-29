% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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

function JR_rot = S5RRRRR2_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_jacobiR_rot_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_jacobiR_rot_4_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:26:51
% EndTime: 2019-03-29 15:26:51
% DurationCPUTime: 0.02s
% Computational Cost: add. (54->16), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->11)
t34 = qJ(1) + qJ(2);
t30 = sin(t34);
t33 = qJ(3) + qJ(4);
t31 = cos(t33);
t36 = t30 * t31;
t29 = sin(t33);
t32 = cos(t34);
t35 = t32 * t29;
t28 = t32 * t31;
t27 = t30 * t29;
t1 = [-t36, -t36, -t35, -t35, 0; t28, t28, -t27, -t27, 0; 0, 0, t31, t31, 0; t27, t27, -t28, -t28, 0; -t35, -t35, -t36, -t36, 0; 0, 0, -t29, -t29, 0; t32, t32, 0, 0, 0; t30, t30, 0, 0, 0; 0, 0, 0, 0, 0;];
JR_rot  = t1;
