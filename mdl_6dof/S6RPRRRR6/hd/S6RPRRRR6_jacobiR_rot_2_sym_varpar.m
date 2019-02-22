% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S6RPRRRR6
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:01
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR6_jacobiR_rot_2_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_jacobiR_rot_2_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_jacobiR_rot_2_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:01:28
% EndTime: 2019-02-22 11:01:28
% DurationCPUTime: 0.02s
% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
t10 = cos(qJ(1));
t9 = sin(qJ(1));
t8 = cos(pkin(11));
t7 = sin(pkin(11));
t1 = [-t9 * t8, 0, 0, 0, 0, 0; t10 * t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t9 * t7, 0, 0, 0, 0, 0; -t10 * t7, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;