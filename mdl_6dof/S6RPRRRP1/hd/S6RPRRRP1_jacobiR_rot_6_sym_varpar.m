% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:51
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:51:04
% EndTime: 2019-02-22 10:51:04
% DurationCPUTime: 0.04s
% Computational Cost: add. (78->18), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->23)
t112 = qJ(3) + qJ(4);
t109 = sin(t112);
t113 = sin(qJ(5));
t120 = t109 * t113;
t114 = cos(qJ(5));
t119 = t109 * t114;
t110 = cos(t112);
t105 = t110 * t113;
t106 = t110 * t114;
t111 = qJ(1) + pkin(10);
t107 = sin(t111);
t118 = t107 * t120;
t117 = t107 * t119;
t108 = cos(t111);
t116 = t108 * t120;
t115 = t108 * t119;
t104 = t108 * t110;
t103 = t107 * t110;
t102 = t108 * t106 + t107 * t113;
t101 = t108 * t105 - t107 * t114;
t100 = t107 * t106 - t108 * t113;
t99 = -t107 * t105 - t108 * t114;
t1 = [-t100, 0, -t115, -t115, -t101, 0; t102, 0, -t117, -t117, t99, 0; 0, 0, t106, t106, -t120, 0; -t107 * t109, 0, t104, t104, 0, 0; t108 * t109, 0, t103, t103, 0, 0; 0, 0, t109, t109, 0, 0; t99, 0, -t116, -t116, t102, 0; t101, 0, -t118, -t118, t100, 0; 0, 0, t105, t105, t119, 0;];
JR_rot  = t1;
