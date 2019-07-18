% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
%
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRRR3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobiR_rot_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobiR_rot_5_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:19:38
% EndTime: 2019-07-18 17:19:38
% DurationCPUTime: 0.11s
% Computational Cost: add. (99->21), mult. (66->20), div. (0->0), fcn. (114->6), ass. (0->27)
t116 = qJ(4) + qJ(5);
t112 = sin(t116);
t117 = qJ(2) + qJ(3);
t113 = sin(t117);
t128 = t113 * t112;
t114 = cos(t116);
t127 = t113 * t114;
t115 = cos(t117);
t126 = t115 * t112;
t118 = sin(qJ(1));
t125 = t118 * t113;
t124 = t118 * t114;
t110 = t118 * t115;
t119 = cos(qJ(1));
t123 = t119 * t113;
t122 = t119 * t114;
t111 = t119 * t115;
t121 = t113 * t124;
t120 = t113 * t122;
t109 = t115 * t114;
t108 = t112 * t123;
t107 = t112 * t125;
t106 = t114 * t111 + t118 * t112;
t105 = -t112 * t111 + t124;
t104 = -t114 * t110 + t119 * t112;
t103 = t112 * t110 + t122;
t1 = [t104, -t120, -t120, t105, t105; t106, -t121, -t121, -t103, -t103; 0, t109, t109, -t128, -t128; t103, t108, t108, -t106, -t106; t105, t107, t107, t104, t104; 0, -t126, -t126, -t127, -t127; -t125, t111, t111, 0, 0; t123, t110, t110, 0, 0; 0, t113, t113, 0, 0;];
JR_rot  = t1;
