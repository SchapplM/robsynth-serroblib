% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:28
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR10_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:28:06
% EndTime: 2019-02-22 11:28:06
% DurationCPUTime: 0.08s
% Computational Cost: add. (55->16), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->26)
t106 = sin(pkin(6));
t108 = sin(qJ(2));
t121 = t106 * t108;
t109 = sin(qJ(1));
t120 = t106 * t109;
t110 = cos(qJ(2));
t119 = t106 * t110;
t111 = cos(qJ(1));
t118 = t106 * t111;
t117 = t109 * t108;
t116 = t109 * t110;
t115 = t111 * t108;
t114 = t111 * t110;
t107 = cos(pkin(6));
t100 = t107 * t115 + t116;
t105 = pkin(11) + qJ(4);
t103 = sin(t105);
t104 = cos(t105);
t113 = t100 * t103 + t104 * t118;
t112 = t100 * t104 - t103 * t118;
t102 = -t107 * t117 + t114;
t101 = t107 * t116 + t115;
t99 = t107 * t114 - t117;
t98 = t102 * t104 + t103 * t120;
t97 = t102 * t103 - t104 * t120;
t1 = [t99, t102, 0, 0, 0, 0; t101, t100, 0, 0, 0, 0; 0, t121, 0, 0, 0, 0; t112, t101 * t104, 0, t97, 0, 0; -t98, -t99 * t104, 0, t113, 0, 0; 0, -t104 * t119, 0, t103 * t121 - t107 * t104, 0, 0; -t113, -t101 * t103, 0, t98, 0, 0; t97, t99 * t103, 0, t112, 0, 0; 0, t103 * t119, 0, t107 * t103 + t104 * t121, 0, 0;];
JR_rot  = t1;
