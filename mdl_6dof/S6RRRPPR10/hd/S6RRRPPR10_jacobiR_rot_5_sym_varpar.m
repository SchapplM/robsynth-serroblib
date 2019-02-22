% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:55
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR10_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:55:38
% EndTime: 2019-02-22 11:55:38
% DurationCPUTime: 0.11s
% Computational Cost: add. (51->24), mult. (163->57), div. (0->0), fcn. (238->10), ass. (0->30)
t107 = sin(pkin(11));
t111 = sin(qJ(3));
t127 = t107 * t111;
t108 = sin(pkin(6));
t126 = t108 * t111;
t114 = cos(qJ(3));
t125 = t108 * t114;
t116 = cos(qJ(1));
t124 = t108 * t116;
t109 = cos(pkin(11));
t123 = t109 * t111;
t115 = cos(qJ(2));
t122 = t111 * t115;
t112 = sin(qJ(2));
t113 = sin(qJ(1));
t121 = t113 * t112;
t120 = t113 * t115;
t119 = t116 * t112;
t118 = t116 * t115;
t110 = cos(pkin(6));
t104 = t110 * t119 + t120;
t99 = -t104 * t111 - t114 * t124;
t117 = -t104 * t114 + t111 * t124;
t106 = -t110 * t121 + t118;
t105 = t110 * t120 + t119;
t103 = t110 * t118 - t121;
t102 = t110 * t111 + t112 * t125;
t101 = t106 * t114 + t113 * t126;
t100 = t106 * t111 - t113 * t125;
t1 = [t103 * t109 + t99 * t107, -t105 * t127 + t106 * t109, t101 * t107, 0, 0, 0; t100 * t107 + t105 * t109, t103 * t127 + t104 * t109, -t117 * t107, 0, 0, 0; 0 (t107 * t122 + t109 * t112) * t108, t102 * t107, 0, 0, 0; -t103 * t107 + t99 * t109, -t105 * t123 - t106 * t107, t101 * t109, 0, 0, 0; t100 * t109 - t105 * t107, t103 * t123 - t104 * t107, -t117 * t109, 0, 0, 0; 0 (-t107 * t112 + t109 * t122) * t108, t102 * t109, 0, 0, 0; t117, -t105 * t114, -t100, 0, 0, 0; t101, t103 * t114, t99, 0, 0, 0; 0, t115 * t125, t110 * t114 - t112 * t126, 0, 0, 0;];
JR_rot  = t1;
