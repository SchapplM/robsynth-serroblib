% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR14V3_jacobiR_rot_2_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobiR_rot_2_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobiR_rot_2_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:04
% EndTime: 2019-04-12 15:12:04
% DurationCPUTime: 0.02s
% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
t7 = sin(qJ(2));
t8 = sin(qJ(1));
t14 = t8 * t7;
t9 = cos(qJ(2));
t13 = t8 * t9;
t10 = cos(qJ(1));
t12 = t10 * t7;
t11 = t10 * t9;
t1 = [-t13, -t12, 0, 0, 0, 0; t11, -t14, 0, 0, 0, 0; 0, t9, 0, 0, 0, 0; t14, -t11, 0, 0, 0, 0; -t12, -t13, 0, 0, 0, 0; 0, -t7, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
