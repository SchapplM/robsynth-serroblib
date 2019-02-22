% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR5
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

function JR_rot = S6RPRRRR5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:00:51
% EndTime: 2019-02-22 11:00:52
% DurationCPUTime: 0.09s
% Computational Cost: add. (135->21), mult. (66->20), div. (0->0), fcn. (114->6), ass. (0->27)
t114 = pkin(11) + qJ(3) + qJ(4);
t112 = sin(t114);
t117 = qJ(5) + qJ(6);
t115 = sin(t117);
t128 = t112 * t115;
t116 = cos(t117);
t127 = t112 * t116;
t113 = cos(t114);
t126 = t113 * t115;
t118 = sin(qJ(1));
t125 = t118 * t115;
t124 = t118 * t116;
t119 = cos(qJ(1));
t123 = t119 * t115;
t122 = t119 * t116;
t121 = t112 * t124;
t120 = t112 * t122;
t111 = t119 * t113;
t110 = t118 * t113;
t109 = t113 * t116;
t108 = t112 * t123;
t107 = t112 * t125;
t106 = t113 * t122 + t125;
t105 = -t113 * t123 + t124;
t104 = -t113 * t124 + t123;
t103 = t113 * t125 + t122;
t1 = [t104, 0, -t120, -t120, t105, t105; t106, 0, -t121, -t121, -t103, -t103; 0, 0, t109, t109, -t128, -t128; t103, 0, t108, t108, -t106, -t106; t105, 0, t107, t107, t104, t104; 0, 0, -t126, -t126, -t127, -t127; -t118 * t112, 0, t111, t111, 0, 0; t119 * t112, 0, t110, t110, 0, 0; 0, 0, t112, t112, 0, 0;];
JR_rot  = t1;
