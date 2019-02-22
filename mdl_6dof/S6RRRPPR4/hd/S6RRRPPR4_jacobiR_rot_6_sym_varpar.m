% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:51
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:51:32
% EndTime: 2019-02-22 11:51:32
% DurationCPUTime: 0.12s
% Computational Cost: add. (118->27), mult. (148->32), div. (0->0), fcn. (221->8), ass. (0->24)
t104 = sin(qJ(6));
t107 = cos(qJ(6));
t103 = qJ(3) + pkin(10);
t101 = sin(t103);
t102 = cos(t103);
t109 = cos(qJ(1));
t106 = sin(qJ(1));
t108 = cos(qJ(2));
t117 = t106 * t108;
t95 = t101 * t117 + t109 * t102;
t96 = -t109 * t101 + t102 * t117;
t115 = t96 * t104 - t95 * t107;
t105 = sin(qJ(2));
t111 = t101 * t104 + t102 * t107;
t119 = t111 * t105;
t116 = t109 * t108;
t114 = t95 * t104 + t96 * t107;
t97 = t101 * t116 - t106 * t102;
t98 = t106 * t101 + t102 * t116;
t113 = t98 * t104 - t97 * t107;
t91 = t97 * t104 + t98 * t107;
t112 = t101 * t107 - t102 * t104;
t93 = t112 * t105;
t1 = [-t114, -t109 * t119, t113, 0, 0, -t113; t91, -t106 * t119, t115, 0, 0, -t115; 0, t111 * t108, -t93, 0, 0, t93; t115, -t109 * t93, t91, 0, 0, -t91; -t113, -t106 * t93, t114, 0, 0, -t114; 0, t112 * t108, t119, 0, 0, -t119; t106 * t105, -t116, 0, 0, 0, 0; -t109 * t105, -t117, 0, 0, 0, 0; 0, -t105, 0, 0, 0, 0;];
JR_rot  = t1;
