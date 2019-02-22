% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:40
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:40:21
% EndTime: 2019-02-22 11:40:21
% DurationCPUTime: 0.04s
% Computational Cost: add. (135->21), mult. (66->20), div. (0->0), fcn. (114->6), ass. (0->27)
t116 = qJ(2) + pkin(11) + qJ(4);
t114 = sin(t116);
t119 = qJ(5) + qJ(6);
t117 = sin(t119);
t130 = t114 * t117;
t118 = cos(t119);
t129 = t114 * t118;
t115 = cos(t116);
t128 = t115 * t117;
t120 = sin(qJ(1));
t127 = t120 * t117;
t126 = t120 * t118;
t121 = cos(qJ(1));
t125 = t121 * t117;
t124 = t121 * t118;
t123 = t114 * t126;
t122 = t114 * t124;
t113 = t121 * t115;
t112 = t120 * t115;
t111 = t115 * t118;
t110 = t114 * t125;
t109 = t114 * t127;
t108 = t115 * t124 + t127;
t107 = -t115 * t125 + t126;
t106 = -t115 * t126 + t125;
t105 = t115 * t127 + t124;
t1 = [t106, -t122, 0, -t122, t107, t107; t108, -t123, 0, -t123, -t105, -t105; 0, t111, 0, t111, -t130, -t130; t105, t110, 0, t110, -t108, -t108; t107, t109, 0, t109, t106, t106; 0, -t128, 0, -t128, -t129, -t129; -t120 * t114, t113, 0, t113, 0, 0; t121 * t114, t112, 0, t112, 0, 0; 0, t114, 0, t114, 0, 0;];
JR_rot  = t1;
