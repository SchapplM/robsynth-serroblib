% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:44
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR9_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:44:41
% EndTime: 2019-02-22 11:44:42
% DurationCPUTime: 0.12s
% Computational Cost: add. (115->18), mult. (117->30), div. (0->0), fcn. (180->8), ass. (0->28)
t110 = sin(pkin(6));
t112 = sin(qJ(2));
t124 = t110 * t112;
t113 = sin(qJ(1));
t123 = t110 * t113;
t114 = cos(qJ(2));
t122 = t110 * t114;
t115 = cos(qJ(1));
t121 = t110 * t115;
t120 = t113 * t112;
t119 = t113 * t114;
t118 = t115 * t112;
t117 = t115 * t114;
t111 = cos(pkin(6));
t103 = t111 * t118 + t119;
t109 = pkin(12) + qJ(4) + qJ(5);
t107 = sin(t109);
t108 = cos(t109);
t97 = -t103 * t108 + t107 * t121;
t116 = t103 * t107 + t108 * t121;
t105 = -t111 * t120 + t117;
t104 = t111 * t119 + t118;
t102 = t111 * t117 - t120;
t101 = -t111 * t107 - t108 * t124;
t100 = -t107 * t124 + t111 * t108;
t99 = t105 * t108 + t107 * t123;
t98 = -t105 * t107 + t108 * t123;
t1 = [t97, -t104 * t108, 0, t98, t98, 0; t99, t102 * t108, 0, -t116, -t116, 0; 0, t108 * t122, 0, t100, t100, 0; t116, t104 * t107, 0, -t99, -t99, 0; t98, -t102 * t107, 0, t97, t97, 0; 0, -t107 * t122, 0, t101, t101, 0; t102, t105, 0, 0, 0, 0; t104, t103, 0, 0, 0, 0; 0, t124, 0, 0, 0, 0;];
JR_rot  = t1;
