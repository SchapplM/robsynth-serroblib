% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:43
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:43:27
% EndTime: 2019-02-22 10:43:27
% DurationCPUTime: 0.09s
% Computational Cost: add. (107->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->25)
t111 = qJ(3) + qJ(4) + pkin(11);
t108 = cos(t111);
t113 = sin(qJ(6));
t121 = t108 * t113;
t112 = qJ(1) + pkin(10);
t109 = sin(t112);
t120 = t109 * t113;
t114 = cos(qJ(6));
t119 = t109 * t114;
t110 = cos(t112);
t118 = t110 * t113;
t117 = t110 * t114;
t107 = sin(t111);
t116 = t107 * t119;
t115 = t107 * t117;
t106 = t108 * t114;
t105 = t110 * t108;
t104 = t109 * t108;
t103 = t107 * t118;
t102 = t107 * t120;
t101 = t108 * t117 + t120;
t100 = -t108 * t118 + t119;
t99 = -t108 * t119 + t118;
t98 = t108 * t120 + t117;
t1 = [t99, 0, -t115, -t115, 0, t100; t101, 0, -t116, -t116, 0, -t98; 0, 0, t106, t106, 0, -t107 * t113; t98, 0, t103, t103, 0, -t101; t100, 0, t102, t102, 0, t99; 0, 0, -t121, -t121, 0, -t107 * t114; -t109 * t107, 0, t105, t105, 0, 0; t110 * t107, 0, t104, t104, 0, 0; 0, 0, t107, t107, 0, 0;];
JR_rot  = t1;
