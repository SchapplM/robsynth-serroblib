% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR2
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
% Datum: 2019-02-22 10:59
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:58:58
% EndTime: 2019-02-22 10:58:58
% DurationCPUTime: 0.09s
% Computational Cost: add. (137->22), mult. (66->20), div. (0->0), fcn. (114->6), ass. (0->24)
t117 = qJ(5) + qJ(6);
t112 = sin(t117);
t118 = qJ(3) + qJ(4);
t113 = sin(t118);
t123 = t113 * t112;
t114 = cos(t117);
t122 = t113 * t114;
t115 = cos(t118);
t121 = t115 * t112;
t109 = t115 * t114;
t116 = qJ(1) + pkin(11);
t110 = sin(t116);
t120 = t110 * t122;
t111 = cos(t116);
t119 = t111 * t122;
t108 = t111 * t115;
t107 = t110 * t115;
t106 = t111 * t123;
t105 = t110 * t123;
t104 = t111 * t109 + t110 * t112;
t103 = t110 * t114 - t111 * t121;
t102 = -t110 * t109 + t111 * t112;
t101 = t110 * t121 + t111 * t114;
t1 = [t102, 0, -t119, -t119, t103, t103; t104, 0, -t120, -t120, -t101, -t101; 0, 0, t109, t109, -t123, -t123; t101, 0, t106, t106, -t104, -t104; t103, 0, t105, t105, t102, t102; 0, 0, -t121, -t121, -t122, -t122; -t110 * t113, 0, t108, t108, 0, 0; t111 * t113, 0, t107, t107, 0, 0; 0, 0, t113, t113, 0, 0;];
JR_rot  = t1;
