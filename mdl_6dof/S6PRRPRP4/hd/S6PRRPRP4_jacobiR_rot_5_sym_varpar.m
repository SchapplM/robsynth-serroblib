% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:46
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRP4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:46:56
% EndTime: 2019-02-22 09:46:56
% DurationCPUTime: 0.07s
% Computational Cost: add. (51->27), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
t110 = sin(pkin(6));
t114 = sin(qJ(3));
t126 = t110 * t114;
t117 = cos(qJ(3));
t125 = t110 * t117;
t112 = cos(pkin(6));
t115 = sin(qJ(2));
t124 = t112 * t115;
t118 = cos(qJ(2));
t123 = t112 * t118;
t113 = sin(qJ(5));
t122 = t113 * t114;
t121 = t113 * t118;
t116 = cos(qJ(5));
t120 = t114 * t116;
t119 = t116 * t118;
t111 = cos(pkin(10));
t109 = sin(pkin(10));
t107 = t112 * t114 + t115 * t125;
t106 = -t112 * t117 + t115 * t126;
t105 = -t109 * t124 + t111 * t118;
t104 = t109 * t123 + t111 * t115;
t103 = t109 * t118 + t111 * t124;
t102 = t109 * t115 - t111 * t123;
t101 = t105 * t117 + t109 * t126;
t100 = t105 * t114 - t109 * t125;
t99 = t103 * t117 - t111 * t126;
t98 = t103 * t114 + t111 * t125;
t1 = [0, -t104 * t122 + t105 * t116, t101 * t113, 0, t100 * t116 - t104 * t113, 0; 0, -t102 * t122 + t103 * t116, t99 * t113, 0, -t102 * t113 + t98 * t116, 0; 0 (t114 * t121 + t115 * t116) * t110, t107 * t113, 0, t106 * t116 + t110 * t121, 0; 0, -t104 * t120 - t105 * t113, t101 * t116, 0, -t100 * t113 - t104 * t116, 0; 0, -t102 * t120 - t103 * t113, t99 * t116, 0, -t102 * t116 - t98 * t113, 0; 0 (-t113 * t115 + t114 * t119) * t110, t107 * t116, 0, -t106 * t113 + t110 * t119, 0; 0, -t104 * t117, -t100, 0, 0, 0; 0, -t102 * t117, -t98, 0, 0, 0; 0, t118 * t125, -t106, 0, 0, 0;];
JR_rot  = t1;
