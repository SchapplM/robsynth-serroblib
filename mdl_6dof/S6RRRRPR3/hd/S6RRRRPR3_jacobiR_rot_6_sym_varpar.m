% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:18
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:18:50
% EndTime: 2019-02-22 12:18:50
% DurationCPUTime: 0.04s
% Computational Cost: add. (94->15), mult. (64->20), div. (0->0), fcn. (111->6), ass. (0->24)
t110 = qJ(2) + qJ(3) + qJ(4);
t108 = sin(t110);
t112 = sin(qJ(1));
t120 = t112 * t108;
t111 = sin(qJ(6));
t119 = t112 * t111;
t113 = cos(qJ(6));
t118 = t112 * t113;
t114 = cos(qJ(1));
t117 = t114 * t108;
t116 = t114 * t111;
t115 = t114 * t113;
t109 = cos(t110);
t107 = t108 * t113;
t106 = t108 * t111;
t105 = t109 * t115;
t104 = t109 * t116;
t103 = t109 * t118;
t102 = t109 * t119;
t101 = -t108 * t119 + t115;
t100 = t108 * t118 + t116;
t99 = t108 * t116 + t118;
t98 = t108 * t115 - t119;
t1 = [t101, t104, t104, t104, 0, t98; t99, t102, t102, t102, 0, t100; 0, t106, t106, t106, 0, -t109 * t113; -t100, t105, t105, t105, 0, -t99; t98, t103, t103, t103, 0, t101; 0, t107, t107, t107, 0, t109 * t111; -t112 * t109, -t117, -t117, -t117, 0, 0; t114 * t109, -t120, -t120, -t120, 0, 0; 0, t109, t109, t109, 0, 0;];
JR_rot  = t1;
